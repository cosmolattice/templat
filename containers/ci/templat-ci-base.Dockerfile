# TempLat CI image -- full rebuild from scratch.
#
# `-DTEMPLAT_TEST=ON` generates one executable per tests/TempLat/**/*.cpp -- 244
# of them, each a C++20 template TU statically linked against Kokkos. A cold
# build of that measures 201 s (serial) / 247 s (mpi) at -j8 on a Zen 4 laptop,
# projecting to ~11 / ~14 min on a 4-vCPU GitHub runner -- not catastrophic, but
# pure waste to repeat on every push to every pull request.
#
# So this image bakes the *build tree*, not just the dependencies: /opt/templat
# holds a git worktree plus a fully populated CMake build directory. A PR job
# checks its merge ref out inside that worktree and runs an incremental `make`,
# which rebuilds only the targets whose include closure actually changed and --
# crucially -- relinks only those. That second part is what a bare ccache cannot
# do: ccache never caches a link step, so every run would still relink all 244
# binaries. Measured here: editing a leaf header rebuilds and relinks exactly 1
# target; editing util/tdd/tdd.h rebuilds all of them.
#
# This file is the expensive path. ci-image.yml normally uses
# templat-ci-update.Dockerfile instead, and only falls back here on the weekly
# cron, on manual request, or when the build configuration itself changed.
#
# Build from the repository root as context:
#   docker buildx build -t ghcr.io/cosmolattice/templat-ci:serial \
#       -f containers/ci/templat-ci-base.Dockerfile \
#       --build-arg config=serial --build-arg threads=8 .

FROM ubuntu:24.04

# serial -> Kokkos Serial backend, no MPI. mpi -> MPI + OpenMP + ParaFaFT.
ARG config=serial
ARG ref=main
ARG threads=4
ARG repo=https://github.com/cosmolattice/templat.git
ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.source=https://github.com/cosmolattice/templat
LABEL org.opencontainers.image.description="TempLat pre-built CI image (Ubuntu 24.04) -- source and populated build tree under /opt/templat"

# Ubuntu 24.04 ships gcc 13 and cmake 3.28; TempLat needs C++20 and (via the
# fetched Kokkos 5.1.1 / KokkosFFT) cmake >= 3.22, so the stock toolchain is
# fine. libfftw3-dev carries the threads and omp variants that
# cmake/libs/fftw/find.cmake looks for; the MPI variants are a separate package.
RUN apt-get -y update \
    && apt-get -y install --no-install-recommends \
        build-essential cmake git ca-certificates ccache libfftw3-dev \
    && if [ "${config}" = "mpi" ]; then \
         apt-get -y install --no-install-recommends \
           libfftw3-mpi-dev openmpi-bin libopenmpi-dev; \
       fi \
    && rm -rf /var/lib/apt/lists/*

# ccache is a second line of defence behind the baked build tree: it catches the
# case where make decides to rebuild a TU whose preprocessed content is actually
# unchanged (a rebase or a `touch` moves mtimes without changing bytes). The
# sloppiness settings are what let it hit across those mtime changes.
ENV CCACHE_DIR=/opt/templat/ccache \
    CCACHE_MAXSIZE=2G \
    CCACHE_COMPILERCHECK=content \
    CCACHE_SLOPPINESS=include_file_mtime,include_file_ctime,time_macros

# OpenMPI refuses to run as root without these, and GitHub Actions container
# jobs *are* root. They have to be environment variables: the obvious
# alternative, threading --allow-run-as-root through MPIEXEC_PREFLAGS, breaks
# because tests/CMakeLists.txt interpolates that variable into an add_test
# `bash -c` string, where a CMake list separator would appear literally as ';'.
# OMP_NUM_THREADS=1 keeps 4 MPI ranks from each spawning a full thread pool.
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMP_NUM_THREADS=1 \
    TEMPLAT_CI_CONFIG=${config}

# A full clone (the repo is ~2.4 MB packed), not a shallow one: PR jobs fetch
# refs/pull/N/merge into this worktree, and git must be able to check that out
# rewriting only the files whose content genuinely differs -- that is precisely
# what keeps mtimes, and therefore make's incremental decisions, correct.
RUN git clone "${repo}" /opt/templat/src \
    && git -C /opt/templat/src checkout --force "${ref}" \
    && git config --global --add safe.directory /opt/templat/src

# The build tree lives *outside* the source tree so that `git clean -xdff` in
# /opt/templat/src -- which PR jobs run to discard any stale generated files --
# cannot destroy it. Both paths are fixed and absolute because CMake build trees
# are not relocatable; PR jobs must build here and not in $GITHUB_WORKSPACE.
#
# NATIVE=OFF is not optional. cmake/system.cmake defaults it ON, which adds
# -march=native in Release. This image is built on one runner and pulled onto
# arbitrary others (GitHub mixes Intel and AMD hosts), so baking the build
# machine's ISA into the objects would SIGILL elsewhere.
RUN if [ "${config}" = "mpi" ]; then \
      extra="-DMPI=ON -DOPENMP=ON -DNPROCESSES=4 -DMPIEXEC_PREFLAGS=--oversubscribe"; \
    else \
      extra="-DNOTHREADING=ON"; \
    fi \
    && cmake -S /opt/templat/src -B /opt/templat/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DTEMPLAT_TEST=ON \
        -DNATIVE=OFF \
        -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
        ${extra} \
    && ccache --zero-stats \
    && cmake --build /opt/templat/build -j "${threads}" \
    && ccache --show-stats

WORKDIR /opt/templat/build
