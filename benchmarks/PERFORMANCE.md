# TempLat performance notes

What is known about where TempLat's time goes on CPU, what has been tried, and what has been
*refuted*. Written down because the ground kept getting re-walked: several of the ideas below are
plausible, were implemented, and lost — and without a record of that, they look attractive again six
months later.

Everything here is measured on the SU(2) pure-gauge leapfrog (`benchmarks/su2_evolution.cpp`), which
is the most arithmetic-heavy kernel in the library and shares a deterministic plane-wave initial
condition with `bench-hila-templat`'s `bench_su2.cpp`, so it doubles as an exact correctness gate.

---

## 1. The lattice kernels must vectorize, and nothing but `test-vectorization` will tell you if they stop

**This is the single largest CPU effect known in TempLat.** Until it was fixed, every lattice kernel
compiled to **100% scalar code** — the SU(2) kick contained 220 scalar FP instructions and *zero*
packed ones. On a machine with two 512-bit FMA units per core that leaves a factor of ~8 unused, and
it is why the kick ran at ~2 GFLOP/s per core, a couple of percent of vector peak.

The cause was not the expression templates and not the Kokkos `View` accessor. It was that **GCC
could not prove the field being written does not alias the fields being read**, so it refused to
vectorize the innermost loop. Kokkos has an escape hatch for exactly this (`KOKKOS_ENABLE_IVDEP_MDRANGE`),
but it only ever defines it for the Intel compiler, so on GCC/Clang its MDRange tile loops get
nothing.

The fix (`parallel/devices/kokkos/`):

- `foreach(name, LayoutStruct, functor)` — the overload every field assignment goes through — now
  dispatches on CPU over the **outer** `NDim-1` dimensions only (`getLocalKokkosOuterPolicy`), and
  walks the contiguous last dimension in a plain loop inside `KokkosNDLambdaWrapperInnerLoop`.
- That inner loop carries `TEMPLAT_ASSUME_INDEPENDENT` (`#pragma GCC ivdep` / `clang loop
  vectorize(assume_safety)`), which asserts what `foreach` already promises its callers: a functor
  dispatched over the lattice is site-independent. Kernels that genuinely accumulate across sites
  must say so with atomics (the radial projector does), and no vectorizer touches those.
- The GPU keeps the full-rank MDRange with the reversed access pattern — that is what gives coalesced
  access there, and it is untouched (`reverse_access_pattern == true` selects it).

The functor is copied **once**, at construction of the wrapper. Do not copy it per row: the
expression trees it holds carry a `shared_ptr` to the memory manager in every leaf, so a copy means
atomic refcount traffic, and copying one per row costs more than vectorizing saves (this was measured
— it was 3x *slower*).

### The measurement

Single core, pinned, best-of-N (this laptop varies ±30% run to run for the *same* binary, so a single
timing run cannot detect even a large regression — always interleave and take the minimum).

| | kick ns/site | drift ns/site | step ns/site |
|---|---|---|---|
| before, 64³ | 125.5 | 37.8 | 163.3 |
| after, 64³ | **46.2** | **9.7** | **55.9**  (2.9x) |
| before, 192³ | 121.6 | 37.8 | 159.4 |
| after, 192³ | **54.5** | **9.8** | **64.3**  (2.5x) |

Packed FP instructions in the CPU lattice kernels: **8 → 8064**. The energy drift is bit-identical
before and after (64³: `plaq = 0.00394059299`, `edrift = 2.696e-03`; 192³: `plaq = 0.02680061607`,
`edrift = 2.245e-03`), as it must be — this changes traversal order, not arithmetic.

### The guardrail

`ctest -R test-vectorization` (`benchmarks/check_vectorization.sh`) disassembles `bench-su2_evolution`
and asserts the CPU lattice kernels contain packed FP instructions. **A scalar kernel passes every
other test in the suite** — it computes exactly the right answer, just ~2.5x too slowly — which is
precisely how this survived as long as it did. It fails on both regressions that matter: the pragma
no longer taking effect, and the CPU dispatch being refactored away.

### How to localise a future regression

Build `benchmarks/kick_ladder.cpp` at `-DKICK_VARIANT=0..3`. It runs the same kick, on the same
layout, at four levels: TempLat as-is (0), a hand-written site loop over `DoEval::eval` (1), the same
with raw-pointer leaves (2), and a hand-written ceiling with no Kokkos and no expression templates (3).
Whichever rung recovers the time tells you which layer is at fault. Variant 1 at 45.8 ns/site vs
variant 0 at 125.5 is what proved the dispatch — not the expression templates — was the problem;
variant 1 with the pragma removed goes straight back to 126.4 ns/site and 0 packed ops, which is what
proved it was the aliasing assertion specifically and not the loop shape.

The ceiling (variant 3) is ~27 ns/site, so there is still ~2x on the table between TempLat and a
hand-written kernel. That gap is *not* the dispatch and not the pragma. Untested candidates: the
`Kokkos::View` accessor's index arithmetic (variant 2 exists to measure this, and has not been run —
it needs a raw-pointer leaf in `ConfigView::eval`), and row alignment (with `nGhosts = 1` every
contiguous row starts at element offset 1, so every vector load is unaligned; padding the contiguous
dimension to a 64 B multiple in `views/fieldviewconfig.h` is cheap to try).

---

## 2. Refuted ideas

Do not re-try these without reading why they lost.

**`Kokkos::MemoryTraits<Unmanaged | Restrict>` on the views.** Zero effect on codegen. The trait does
not reach the vectorizer as a `__restrict__` on the data pointer.

**`-fvect-cost-model=unlimited`.** Still scalar, and 7x slower overall. GCC was not *declining* to
vectorize on cost grounds; it *could not*, because of the aliasing it could not disprove. Fixing the
aliasing is what worked.

**The SU(2) scatter kick** (`-DKICK_SCATTER` in `su2_evolution.cpp`; the derivation is in that file's
header). Each site visits 3 unique plaquettes instead of 12, cutting arithmetic ~4x, and it is
energy-exact. It is nonetheless **3.2–4x slower** than the gather kick (201.9 vs 50.1 ns/site at 64³;
199.5 vs 61.7 at 192³) and it will stay that way: it scatters read-modify-writes to Π(x) *and*
Π(x+ĵ), so along the contiguous dimension it writes both Π(x) and Π(x+ẑ) — a genuine loop-carried
dependence. It therefore cannot take the `ivdep` path above (asserting independence there would be
*incorrect*, not merely optimistic), and it stays scalar while the gather kick vectorizes. This was
already 1.4–1.7x slower against the *scalar* gather, for a different reason (L1-bound on the
scattered RMW, because the SU(2) components live in separate `Field`s). Fewer FLOPs is the wrong
thing to optimise for here.

**Hoisting the contiguous dimension with a per-row stack copy of the functor.** 3x slower. See the
note about `shared_ptr` leaves above. The current fix hoists the dimension *without* the copy.

---

## 3. Where the time goes now

- The **kick is compute-bound** (~98% of its time is arithmetic, not memory traffic) and is ~85% of a
  leapfrog step. It was the right thing to attack, and vectorization was the right lever.
- The **drift is bandwidth-bound**, and is now ~15% of the step. Its exp map (`su2algebra/helpers/
  paulivectorsalgebra.h`) contains an `a > 1e-15` sinc guard — a data-dependent branch in the loop
  body, which normally blocks vectorization outright. GCC if-converts it here (the drift vectorized:
  37.8 → 9.7 ns/site), so it is not currently costing anything, but it is fragile: if the drift ever
  goes scalar again, replace the branch with a branchless select before looking anywhere else.
- **Kick/drift sweep fusion** (removing the drift's DRAM round-trip of Π) was refuted against the
  scalar baseline. It is worth strictly less now: the drift is a smaller share of the step than it
  was, so the ceiling on this idea is ~10%. It is also not free to make correct — a fused sweep
  updates U(x) while neighbouring sites still need the *old* U(x) for their kick.
- **Hybrid MPI + threading loses**, and this is unaffected by the vectorization work: every
  expression-template assignment is its own fork/join, so the thread team is re-formed per kernel.
- The FFT forces a 2D pencil decomposition, which constrains the domain decomposition independently
  of any of the above.

---

## Reproducing

```bash
cmake -B build -DTEMPLAT_TEST=ON -DTEMPLAT_BENCH=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8 --target bench-su2_evolution
ctest --test-dir build -R test-vectorization        # codegen guardrail

# ns/site, pinned, and the correctness gate (plaq / edrift must not move)
taskset -c 2 ./build/benchmarks/bench-su2_evolution

# why is it slow again? -- the ladder
c++ -O3 -march=native -DKICK_VARIANT=1 benchmarks/kick_ladder.cpp ...
# and to see what the vectorizer refused, and where:
c++ -O3 -march=native -fopt-info-vec-missed ...
```

Cap builds at `-j 8` on this machine.
