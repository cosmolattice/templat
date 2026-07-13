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

Build `benchmarks/kick_ladder.cpp` at `-DKICK_VARIANT=0..5`. It runs the same kick, on the same layout,
at rungs of increasing distance from the library: TempLat as-is (0), a hand-written site loop over
`DoEval::eval` (1), raw-pointer leaves (2 — refuted, see below), a hand-written ceiling with no Kokkos
and no expression templates (3), and the de-fused kick (4, 5 — see §2).
Whichever rung recovers the time tells you which layer is at fault. Variant 1 at 45.8 ns/site vs
variant 0 at 125.5 is what proved the dispatch — not the expression templates — was the problem;
variant 1 with the pragma removed goes straight back to 126.4 ns/site and 0 packed ops, which is what
proved it was the aliasing assertion specifically and not the loop shape.

---

## 2. The remaining gap to a hand-written kernel is register spills

The ceiling (variant 3) is ~27–32 ns/site against TempLat's ~46–61, so there is still ~1.5–2x between
the library and a hand-written kernel. It has been chased, and it is **not** the dispatch, **not** the
pragma, and **not** the leaf accessor. It is that the fused kick expression does not fit in the
register file.

`perf stat`, 64³, single core (TempLat runs kick+drift, the ceiling only the kick, so the FP counts
are the fair comparison, not the totals):

| | instructions | cycles | FP ops | IPC | loads | stores |
|---|---|---|---|---|---|---|
| TempLat | 10.7e9 | 7.77e9 | 33.0e9 | 1.38 | 5.33e9 | 1.61e9 |
| ceiling | 3.64e9 | 3.22e9 | 26.0e9 | 1.13 | 1.22e9 | 0.40e9 |

TempLat's **IPC is higher** than the ceiling's, and it does comparable FP work. It simply executes ~3x
the instructions to do it: **4.4x the loads and 4x the stores.** Statically, the kick loop contains
**3276 stack-touching vector moves against the ceiling's 83.** `Pi(i) = Pi(i) - dt * Total(j, plaq -
plaqBack)` materializes all twelve four-link products in one expression; vectorized, each live
quaternion component is a full zmm, and the live set does not fit in 32 registers. GCC spills it and
reloads. (It is not unrolling: `-fno-unroll-loops` leaves all 3276 spills in place.)

**The lever is to stop fusing so much.** Splitting the kick into one assignment per transverse
direction — `for j != i: Pi(i) = Pi(i) - dt*(plaq(U,i,j) - plaqBack(U,i,j))` — trades one extra Π
read+write for a live set that fits, and is **11–22% faster**, with `plaq`/`edrift` bit-identical.
Spills drop 3276 -> 1500. `su2_evolution.cpp` **now does this by default**; `-DKICK_FUSED` restores the
fused form for comparison.

| kick, ns/site | fused | de-fused |
|---|---|---|
| 64³ | 46.9 | **41.9** |
| 192³ | 55.9 | **43.4** |

Do not over-split: one assignment per *term* (`kick_ladder` variant 5) drops spills further, to 1166,
but is slower — the extra Π traffic then outweighs the spills it saves. There is an optimum and it is
not at either end.

**This applies to CosmoLattice, and probably more so.** `CosmoInterface/evolvers/kernels/su2kernels.h`
builds `-normGrad * (GradSU2 - GradSU2Back) - normSU2Source * SU2Source`, where `GradSU2` and
`GradSU2Back` are each a `Total(j, ...)` over plaquettes, and the evolvers apply it as a single
`piSU2(n) += (w * dt) * SU2Kernels::get(model, n)`. That is the fused shape with an *extra* source term
in the live set. Not yet measured there — but it is the same expression pathology, and the same
de-fusing should apply.

Closing the rest of the gap means shrinking the live set *inside* the plaquette evaluation — a fused
primitive that accumulates the algebra components without materializing every intermediate quaternion
(compare the `su2dotter` fused-eval primitive, which was introduced for the same class of reason).
That is real work and it has not been attempted.

### Two things that sound right and are not

**Raw-pointer leaves (`RawAccessor`) — built, measured, refuted.** Replacing `Kokkos::View::operator()`
in `ConfigView::eval` with a POD `T* __restrict__` + strides accessor performs **identically to the
View** (52.6 vs 53.8 ns/site at 64³; 61.9 vs 59.7 at 192³) and generates near-identical code. The View
accessor is not the cost. The code was reverted rather than left behind a flag.

If you do rebuild it, note the trap that cost the first attempt: the **last stride must be a
compile-time 1**, not a stored runtime value. Storing it costs *3x the packed FP instructions and 2x
the integer multiplies*, because unit-stride along the contiguous dimension is exactly what tells the
vectorizer a row can be loaded as a vector instead of gathered. `Kokkos::View` gets this right for
`LayoutRight`; a naive replacement does not.

**Row alignment.** With `nGhosts = 1` every contiguous row starts at element offset 1, so every vector
load is unaligned — but the hand-written ceiling has exactly the same property, so this *cannot*
explain the gap between them. It is a candidate for beating the ceiling, not for reaching it.

**A kernel-form snapshot of the tree (`asKernel()`, POD leaves).** This was the planned fallback for
this gap. It would make leaf pointers loop-invariant — which the measurements above say was never the
problem. It does nothing about register pressure. Do not spend the ~30 operator node types on it.

---

## 3. Refuted ideas

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

## 4. Where the time goes now

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
