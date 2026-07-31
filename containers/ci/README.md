# CI container images

`-DTEMPLAT_TEST=ON` generates **one executable per `tests/TempLat/**/*.cpp`** — 244 of
them, each a C++20 template TU statically linked against Kokkos. A cold build of all of
that measures **201 s** (serial) and **247 s** (mpi) at `-j8` on a Zen 4 laptop, which
projects to roughly **11 min** and **14 min** on a 4-vCPU GitHub runner. Not catastrophic,
but paying it on every push to every pull request is pure waste.

These images make it unnecessary. Each one carries, under `/opt/templat`:

| path | what it is |
|---|---|
| `/opt/templat/src` | a full git worktree of TempLat, `origin` pointing at GitHub |
| `/opt/templat/build` | a fully configured **and built** CMake tree — all 244 test binaries |
| `/opt/templat/ccache` | a warm ccache |

A PR job moves the worktree to its merge ref and runs an incremental `make`. Only the targets
whose include closure actually changed get rebuilt, and — the part a bare ccache cannot give
you — only those get **relinked**.

Measured on the serial image, over a 61-target slice, `-j8`:

| what changed | `make` | compiled | relinked |
|---|---|---|---|
| nothing | 5 s | 0 | 0 |
| `touch` only, content identical | 6 s | 1 | 1 |
| a leaf header (`util/spline.h`) | 6 s | 1 | 1 |
| a core header (`util/tdd/tdd.h`) | 34 s | 122 | 61 |

The `touch` row is why the ccache is in the image: a rebase or a mtime-only change makes
`make` re-run the compiler, and ccache turns that back into a no-op.

Running the suite is cheap by comparison — all 244 tests on 4 CPUs: **19 s** serial
(`ctest -j4`), **125 s** mpi (`ctest -j1`, 4 ranks per test, 1.2 GiB peak RSS).

Two configs are published:

| tag | flags |
|---|---|
| `ghcr.io/cosmolattice/templat-ci:serial` | `-DNOTHREADING=ON` |
| `ghcr.io/cosmolattice/templat-ci:mpi` | `-DMPI=ON -DOPENMP=ON -DNPROCESSES=4 -DMPIEXEC_PREFLAGS=--oversubscribe` |

both with `-DCMAKE_BUILD_TYPE=Release -DTEMPLAT_TEST=ON -DNATIVE=OFF
-DCMAKE_CXX_COMPILER_LAUNCHER=ccache`. `MPI=ON` force-enables `PARAFAFT`, so the mpi image
also bakes the ParaFaFT clone — PR jobs are then immune to a GitHub outage at configure time.

## How they are refreshed

`.github/workflows/ci-image.yml`, the one workflow that does *not* run on pull requests:

- **push to `main`** → incremental (`templat-ci-update.Dockerfile`), a few minutes.
- **weekly cron, `workflow_dispatch` with `full=true`, or a change under `cmake/`,
  `CMakeLists.txt`, `tests/CMakeLists.txt`, `containers/ci/`** → full rebuild
  (`templat-ci-base.Dockerfile`), ~15–25 min including apt and the Kokkos/ParaFaFT fetch.

If an incremental refresh fails, the workflow falls back to a full rebuild in the same run,
so a broken increment cannot leave a stale tag silently slowing every PR down.

**The weekly full rebuild is load-bearing, not hygiene.** Every incremental refresh appends a
Docker layer holding whatever was rebuilt, and a merge touching a core header can add a few
hundred MB. The cron is what bounds that growth to about seven layers. (GHCR storage is free
for public packages, so the cost of growth is pull time, not money.)

## One-time setup

After the first successful run, set the GHCR package **`cosmolattice/templat-ci` to Public**
(Packages → templat-ci → Package settings → Change visibility). `ci.yml` deliberately has no
registry-login step, so that pull requests from forks work; that only holds if the package is
public. A `401` when pulling the image is this setting having reverted.

## Three things that will bite you if changed

1. **`-DNATIVE=OFF` is mandatory.** `cmake/system.cmake` defaults `NATIVE=ON`, adding
   `-march=native` in Release. These images are built on one runner and pulled onto arbitrary
   others — GitHub mixes Intel and AMD hosts — so baking the build machine's ISA into the
   object files would SIGILL elsewhere.
2. **Do not pass `--allow-run-as-root` through `MPIEXEC_PREFLAGS`.** GitHub Actions container
   jobs run as root and OpenMPI refuses to. But `tests/CMakeLists.txt` interpolates
   `MPIEXEC_PREFLAGS` into an `add_test` `bash -c` string, where a CMake list separator would
   appear literally as `;` and break the command. The images set
   `OMPI_ALLOW_RUN_AS_ROOT=1` / `OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` instead. `MPIEXEC_PREFLAGS`
   carries only `--oversubscribe`, a single token — a 4-vCPU runner is typically 2 physical
   cores, so OpenMPI sees 2 slots and would otherwise refuse 4 ranks.
3. **The build tree lives outside the source tree, at fixed absolute paths.** CMake build
   trees are not relocatable, so PR jobs must build at `/opt/templat/build` and *not* in
   `$GITHUB_WORKSPACE` (which a container job mounts at `/__w/templat/templat`). Keeping
   `build/` out of `src/` is also what lets a PR job run `git clean -xdff` in the worktree
   without destroying it.

Related: `cmake/device.cmake` auto-detects the backend when no device flag is given, so both
configs pass one explicitly.

## Working with the images locally

Build one (from the repository root — the Dockerfile clones from GitHub rather than using the
context, so it builds whatever `ref` names, not your working tree):

```bash
docker build -t templat-ci:serial \
    -f containers/ci/templat-ci-base.Dockerfile \
    --build-arg config=serial --build-arg threads=8 --build-arg ref=main .
```

Reproduce what a PR job does, against a branch you have pushed:

```bash
docker run --rm -it --shm-size=2g templat-ci:serial bash -c '
  cd /opt/templat/src
  git fetch --no-tags --force origin "+my-branch:refs/ci/target"
  git checkout --force refs/ci/target && git clean -xdff
  cmake -S /opt/templat/src -B /opt/templat/build
  cd /opt/templat/build
  ctest -N -I 1,,2 | sed -n "s/^ *Test *#[0-9]\+: //p" > /tmp/t.txt
  time make -j8 $(tr "\n" " " < /tmp/t.txt)
  ctest -I 1,,2 --output-on-failure -j8'
```

Check the MPI plumbing (root, oversubscribe, `/dev/shm`) without a full run:

```bash
docker run --rm --shm-size=2g templat-ci:mpi \
    ctest -R test-mpiallreduce --output-on-failure
```

To iterate on uncommitted changes, bind-mount over the worktree instead of cloning — note
that this rewrites every mtime, so the first build inside will be slow unless ccache saves it:

```bash
docker run --rm -it --shm-size=2g -v "$PWD:/opt/templat/src:ro" templat-ci:serial bash
```
