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
  // hot kernel. On CPU we dispatch only over the outer dimensions and walk the contiguous one in the
  // functor, so that the vectorizer gets a loop it can actually reason about; see
  // KokkosNDLambdaWrapperInnerLoop. On GPU the full-rank MDRange is what gives coalesced access, and
  // is kept.
  template <size_t NDim, typename Functor>
    requires requires(Functor functor) { functor(device_kokkos::IdxArray<NDim>{}); }
  void foreach (const std::string &name, const LayoutStruct<NDim> &mLayout, const Functor &functor)
  {
    // -DTEMPLAT_LEGACY_MDRANGE restores the pre-vectorization dispatch (Kokkos' full-rank tiled
    // MDRange, which emits scalar code). It exists to A/B the dispatch in situ: the CPU kernels are
    // 2.5x faster this way on every machine and geometry measured HERE, but a cluster run reported the
    // SU(2) kick 10% SLOWER at scale, and that has not been reproduced or explained. Until it is, keep
    // the escape hatch -- and see benchmarks/PERFORMANCE.md before deleting it.
#ifndef TEMPLAT_LEGACY_MDRANGE
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
