#ifndef TEMPLAT_PARALLEL_KOKKOS_ITERATION_H
#define TEMPLAT_PARALLEL_KOKKOS_ITERATION_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include "TempLat/parallel/devices/kokkos/kokkos.h"
#include "TempLat/lattice/memory/memorylayouts/checkerboardlayout.h"

namespace TempLat::device_kokkos::iteration
{
  // ===== FOREACH =====
  // 1. foreach: device_kokkos::array
  template <size_t NDim, typename Functor, typename I = device_kokkos::Idx>
    requires requires(Functor functor) { functor(device_kokkos::IdxArray<NDim>{}); }
  void foreach (const std::string &name, const device_kokkos::array<I, NDim> &starts,
                const device_kokkos::array<I, NDim> &stops, const Functor &functor)
  {
    Kokkos::parallel_for(name, getLocalKokkosPolicy(starts, stops),
                         device_kokkos::KokkosNDLambdaWrapper<NDim, Functor>(functor));
  }
  // 2. foreach: LayoutStruct. This is the overload every field assignment goes through, i.e. every
  // hot kernel. The default is Kokkos' full-rank tiled MDRange. -DTEMPLAT_VECTORIZED_MDRANGE opts
  // into dispatching over the outer dimensions only and walking the contiguous one inside the
  // functor, which is what lets the vectorizer see the innermost loop (KokkosNDLambdaWrapperInnerLoop).
  //
  // ===========================================================================================
  // WHY THE VECTORIZED DISPATCH IS *NOT* THE DEFAULT -- REVISIT THIS, IT IS NOT SETTLED
  // ===========================================================================================
  // It used to be. It was flipped back after an A/B of the two dispatches (one binary per macro,
  // same source, interleaved reps) reproduced on PC2 Noctua 2 (2x EPYC 7763, Slurm job 33674207)
  // the "+10% kick at scale" that benchmarks/PERFORMANCE.md had recorded as unexplained. Pure-gauge
  // SU(2) leapfrog. Measured kick, ns/site/step, all min-max envelopes disjoint over 3 reps:
  //
  //     N     ranks   vectorized   legacy     winner
  //     384    16        9.82      11.89      vectorized +21%
  //     384    32        9.23       8.67      legacy      +6.4%
  //     384    64        4.22       4.03      legacy      +4.7%
  //     384   128        2.25       2.13      legacy      +5.5%
  //     768   128        2.36       2.26      legacy      +4.6%
  //     768   256        1.20       1.14      legacy      +5.3%
  //
  // The controlling variable is NODE OCCUPANCY, not lattice size and not per-rank volume. The
  // pencil decomposition happens to give two exactly-matched local boxes, and they disagree:
  //   384^3 @ 16 ranks and 768^3 @ 128 ranks are both 96x96x384 = 3.54M sites/rank -- vectorized
  //   wins the first by 21% and loses the second by 4.6%. Ditto 384^3 @ 32 vs 768^3 @ 256 (both
  //   1.77M sites/rank), where legacy wins both. What differs is how many cores are competing for
  // memory bandwidth. Below ~25% occupancy the kernel is compute-bound and removing instructions
  // pays; at full occupancy it is DRAM-bound and only the dispatch's residual cost is visible.
  // Since production runs fill nodes, legacy is the right default and the win above is unreachable.
  //
  // Two facts worth keeping when this is picked up again:
  //   - The penalty is confined to STENCIL kernels. Drift (U = exp(dt*Pi)*U, site-local) is a dead
  //     heat at every occupancy; the kick (the only neighbour-reading kernel) carries all of it.
  //   - It is not "the inner sweep is too long for cache". At fixed 128 ranks, doubling the
  //     contiguous sweep 384 -> 768 does not worsen it (5.5% -> 4.6%). Consistent with the
  //     strip-mining experiment already refuted in PERFORMANCE.md.
  //
  // THE INTENDED FIX is to choose per launch instead of per build, since neither path dominates.
  // The signal is already computed and needs no new dependency:
  //     occupancy = ThreadSettings::getMPILocalSize() * ThreadSettings::getMaxThreadCount()
  //               / <cores on node>            // parallel/threadsettings.h, filled at MPI_Init
  //                                            // from MPI_Comm_split_type(MPI_COMM_TYPE_SHARED)
  // It was not done now because three things need resolving first, and none is a detail:
  //   1. COMPILE TIME. A runtime branch instantiates BOTH wrappers for every kernel in every TU.
  //      An SU2Axion TU is already ~162 s, 64% of it backend inlining -- exactly what this lands
  //      on. Measure before committing. Mitigation: the drift result above says only
  //      neighbour-reading expressions need both paths, so a GhostsHunter-style trait could keep
  //      site-local kernels single-instantiation.
  //   2. THE THRESHOLD IS UNCALIBRATED. The data brackets the crossover only between 12.5% and 25%
  //      occupancy, on ONE machine. The mechanism is bandwidth saturation, so it moves with
  //      bytes-per-core; Otus (2x EPYC 9655, AVX-512) will sit elsewhere and has not been measured.
  //      A hardcoded constant that silently picks wrong is worse than today's explicit flag.
  //   3. hardware_concurrency() IS SMT- AND AFFINITY-BLIND. It reports online logical CPUs, so a
  //      full node reads 0.5 rather than 1.0 when SMT is on, and it ignores cgroup limits.
  // Whatever lands should keep an explicit override, so a wrong heuristic is never silent.
  //
  // GPU is unaffected by the default either way: the opt-in branch is guarded on
  // !reverse_access_pattern, and the full-rank MDRange is what gives coalesced access there.
  // ===========================================================================================
#if defined(TEMPLAT_VECTORIZED_MDRANGE) && defined(TEMPLAT_LEGACY_MDRANGE)
#error "TEMPLAT_VECTORIZED_MDRANGE and TEMPLAT_LEGACY_MDRANGE are mutually exclusive. Legacy is now \
the default, so -DTEMPLAT_LEGACY_MDRANGE is redundant; drop it, or drop the vectorized one."
#endif
  template <size_t NDim, typename Functor>
    requires requires(Functor functor) { functor(device_kokkos::IdxArray<NDim>{}); }
  void foreach (const std::string &name, const LayoutStruct<NDim> &mLayout, const Functor &functor)
  {
#ifdef TEMPLAT_VECTORIZED_MDRANGE
    if constexpr (!device_kokkos::reverse_access_pattern && NDim >= 2) {
      const auto localSizes = mLayout.getSizesInMemory();
      const device_kokkos::Idx nGhosts = mLayout.getNGhosts();
      Kokkos::parallel_for(name, device_kokkos::getLocalKokkosOuterPolicy(mLayout),
                           device_kokkos::KokkosNDLambdaWrapperInnerLoop<NDim, Functor>(
                               functor, nGhosts, nGhosts + (device_kokkos::Idx)localSizes[NDim - 1]));
    } else
#endif
    {
      Kokkos::parallel_for(name, device_kokkos::getLocalKokkosPolicy(mLayout),
                           device_kokkos::KokkosNDLambdaWrapper<NDim, Functor>(functor));
    }
  }
  // 3. foreach: CheckerboardLayout
  template <size_t NDim, typename Functor>
    requires requires(Functor functor) { functor(device_kokkos::IdxArray<NDim>{}); }
  void foreach (const std::string &name, const CheckerboardLayout<NDim> &cb, const Functor &functor)
  {
    CheckerboardForEachWrapper<NDim, Functor> wrapped{cb, functor};
    Kokkos::parallel_for(
        name, getLocalKokkosPolicy(cb.getStarts(), cb.getStops()),
        device_kokkos::KokkosNDLambdaWrapper<NDim, CheckerboardForEachWrapper<NDim, Functor>>(wrapped));
  }

  // ===== REDUCE =====
  // 1. reduce: device_kokkos::array -> value
  template <size_t NDim, typename Functor, typename T, typename I = device_kokkos::Idx>
    requires requires(Functor functor, T &update) { functor(device_kokkos::IdxArray<NDim>{}, update); }
  void reduce(const std::string &name, const device_kokkos::array<I, NDim> &starts,
              const device_kokkos::array<I, NDim> &stops, const Functor &functor, T &result)
  {
    Kokkos::parallel_reduce(name, getLocalKokkosPolicy(starts, stops),
                            device_kokkos::KokkosNDLambdaWrapperReduction<NDim, Functor>(functor), result);
  }
  // 2. reduce: LayoutStruct -> value
  template <size_t NDim, typename Functor, typename T>
    requires requires(Functor functor, T &update) { functor(device_kokkos::IdxArray<NDim>{}, update); }
  void reduce(const std::string &name, const LayoutStruct<NDim> &mLayout, const Functor &functor, T &result)
  {
    Kokkos::parallel_reduce(name, device_kokkos::getLocalKokkosPolicy(mLayout),
                            device_kokkos::KokkosNDLambdaWrapperReduction<NDim, Functor>(functor), result);
  }
  // 4. reduce: device_kokkos::array -> View or Reduction
  template <size_t NDim, typename Functor, typename View, typename I = device_kokkos::Idx>
    requires requires(Functor functor, typename View::value_type &update) {
      functor(device_kokkos::IdxArray<NDim>{}, update);
    }
  void reduce(const std::string &name, const device_kokkos::array<I, NDim> &starts,
              const device_kokkos::array<I, NDim> &stops, const Functor &functor, View view)
  {
    Kokkos::parallel_reduce(name, getLocalKokkosPolicy(starts, stops),
                            device_kokkos::KokkosNDLambdaWrapperReduction<NDim, Functor>(functor), view);
  }
  // 5. reduce: LayoutStruct -> View or Reduction
  template <size_t NDim, typename Functor, typename View>
    requires requires(Functor functor, typename View::value_type &update) {
      functor(device_kokkos::IdxArray<NDim>{}, update);
    }
  void reduce(const std::string &name, const LayoutStruct<NDim> &mLayout, const Functor &functor, View view)
  {
    Kokkos::parallel_reduce(name, device_kokkos::getLocalKokkosPolicy(mLayout),
                            device_kokkos::KokkosNDLambdaWrapperReduction<NDim, Functor>(functor), view);
  }
  // 6. reduce: CheckerboardLayout -> value
  template <size_t NDim, typename Functor, typename T>
    requires requires(Functor functor, T &update) { functor(device_kokkos::IdxArray<NDim>{}, update); }
  void reduce(const std::string &name, const CheckerboardLayout<NDim> &cb, const Functor &functor, T &result)
  {
    CheckerboardReduceWrapper<NDim, Functor> wrapped{cb, functor};
    Kokkos::parallel_reduce(
        name, getLocalKokkosPolicy(cb.getStarts(), cb.getStops()),
        device_kokkos::KokkosNDLambdaWrapperReduction<NDim, CheckerboardReduceWrapper<NDim, Functor>>(wrapped), result);
  }
  // 7. reduce: CheckerboardLayout -> View or Reduction
  template <size_t NDim, typename Functor, typename View>
    requires requires(Functor functor, typename View::value_type &update) {
      functor(device_kokkos::IdxArray<NDim>{}, update);
    }
  void reduce(const std::string &name, const CheckerboardLayout<NDim> &cb, const Functor &functor, View view)
  {
    CheckerboardReduceWrapper<NDim, Functor> wrapped{cb, functor};
    Kokkos::parallel_reduce(
        name, getLocalKokkosPolicy(cb.getStarts(), cb.getStops()),
        device_kokkos::KokkosNDLambdaWrapperReduction<NDim, CheckerboardReduceWrapper<NDim, Functor>>(wrapped), view);
  }

  // ===== REDUCERS =====
  using Kokkos::Max;
  using Kokkos::Min;
  using Kokkos::Prod;
  using Kokkos::Sum;

  // ===== FENCE =====
  inline void fence() { Kokkos::fence(); }
} // namespace TempLat::device_kokkos::iteration

#endif
