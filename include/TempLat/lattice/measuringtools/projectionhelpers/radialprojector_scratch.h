#ifndef TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_RADIALPROJECTOR_SCRATCH_H
#define TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_RADIALPROJECTOR_SCRATCH_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#ifdef DEVICE_KOKKOS

#include "TempLat/parallel/devices/kokkos/kokkos.h"

namespace TempLat
{
  /**
   * @brief Decompose a flat site index into N-dimensional ghost-padded memory indices.
   *
   * Uses row-major (LayoutRight) strides: the last dimension varies fastest.
   * This naturally gives coalesced GPU memory access for LayoutRight views.
   *
   * @param flatIdx  Flat index in [0, product(memSizes))
   * @param memSizes Local sizes in memory order (from LayoutStruct::getSizesInMemory())
   * @param nGhosts  Ghost padding offset (from LayoutStruct::getNGhosts())
   * @return IdxArray<NDim> with each element in [nGhosts, nGhosts + memSizes[d])
   */
  template <size_t NDim>
  DEVICE_FORCEINLINE_FUNCTION
  device_kokkos::IdxArray<NDim> flatToMemoryIndex(
      ptrdiff_t flatIdx,
      const device_kokkos::IdxArray<NDim>& memSizes,
      device_kokkos::Idx nGhosts)
  {
    device_kokkos::IdxArray<NDim> idx;
    for (int d = static_cast<int>(NDim) - 1; d >= 0; --d) {
      idx[d] = nGhosts + (flatIdx % memSizes[d]);
      flatIdx /= memSizes[d];
    }
    return idx;
  }

  /**
   * @brief Manages the layout of 9 bin arrays in Kokkos scratch memory.
   *
   * Arrays (all of type T, size nBins each):
   *   0: values averages      1: values variances
   *   2: values mins           3: values maxs
   *   4: binBounds averages    5: binBounds variances
   *   6: binBounds mins        7: binBounds maxs
   *   8: multiplicities
   *
   * On GPU Level 0: these live in shared memory (~48 KB per block).
   * On GPU Level 1 or CPU: these live in global/heap memory per team.
   */
  template <typename T>
  struct ScratchBinLayout
  {
    static constexpr int kNumArrays = 9;

    using ScratchView = Kokkos::View<T*,
        typename device_kokkos::DefaultExecutionSpace::scratch_memory_space,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    /// Total scratch bytes needed for nBins bins.
    static size_t bytesNeeded(ptrdiff_t nBins)
    {
      return static_cast<size_t>(kNumArrays) * ScratchView::shmem_size(static_cast<size_t>(nBins));
    }

    ptrdiff_t mNBins;

    // Raw pointers into scratch memory for each array.
    T* vAvg;   // values averages
    T* vVar;   // values variances
    T* vMin;   // values mins
    T* vMax;   // values maxs
    T* bAvg;   // binBounds averages
    T* bVar;   // binBounds variances
    T* bMin;   // binBounds mins
    T* bMax;   // binBounds maxs
    T* mult;   // multiplicities

    /**
     * @brief Construct from a Kokkos team scratch memory handle.
     *
     * Uses ScratchView for proper alignment, then stores raw T* pointers
     * for efficient per-element access.
     */
    template <typename ScratchSpace>
    DEVICE_FUNCTION
    ScratchBinLayout(const ScratchSpace& scratch, ptrdiff_t nBins)
        : mNBins(nBins)
    {
      vAvg = ScratchView(scratch, nBins).data();
      vVar = ScratchView(scratch, nBins).data();
      vMin = ScratchView(scratch, nBins).data();
      vMax = ScratchView(scratch, nBins).data();
      bAvg = ScratchView(scratch, nBins).data();
      bVar = ScratchView(scratch, nBins).data();
      bMin = ScratchView(scratch, nBins).data();
      bMax = ScratchView(scratch, nBins).data();
      mult = ScratchView(scratch, nBins).data();
    }

    /// Initialize bin `b` to identity values. Call from TeamThreadRange over nBins.
    DEVICE_FORCEINLINE_FUNCTION
    void init(ptrdiff_t b) const
    {
      vAvg[b] = T(0);
      vVar[b] = T(0);
      vMin[b] = Kokkos::reduction_identity<T>::min();   // +inf (or +max)
      vMax[b] = Kokkos::reduction_identity<T>::max();   // -inf (or -max)
      bAvg[b] = T(0);
      bVar[b] = T(0);
      bMin[b] = Kokkos::reduction_identity<T>::min();
      bMax[b] = Kokkos::reduction_identity<T>::max();
      mult[b] = T(0);
    }

    /// Accumulate one lattice site's contribution into scratch bin `b`.
    /// Uses atomics (fast shared-memory atomics on GPU Level 0; plain writes on CPU with team_size=1).
    DEVICE_FORCEINLINE_FUNCTION
    void accumulate(ptrdiff_t b, const T& value, const T& position, const T& weight) const
    {
      // Values
      Kokkos::atomic_add(&vAvg[b], weight * value);
      Kokkos::atomic_add(&vVar[b], weight * value * value);
      Kokkos::atomic_min(&vMin[b], value);
      Kokkos::atomic_max(&vMax[b], value);
      // Bin bounds
      Kokkos::atomic_add(&bAvg[b], weight * position);
      Kokkos::atomic_add(&bVar[b], weight * position * position);
      Kokkos::atomic_min(&bMin[b], position);
      Kokkos::atomic_max(&bMax[b], position);
      // Multiplicities
      Kokkos::atomic_add(&mult[b], weight);
    }

    /**
     * @brief Merge scratch bin `b` into global device views.
     *
     * Call from TeamThreadRange over nBins after a team_barrier.
     * Uses global-memory atomics (one per bin per team — minimal contention).
     *
     * Skips bins that received no contributions (mult[b] == 0) to avoid
     * polluting min/max with identity values.
     */
    template <typename ViewT, typename ViewM>
    DEVICE_FORCEINLINE_FUNCTION
    void mergeTo(
        // RadialProjectionSingleQuantity device views (values)
        const ViewT& gVAvg, const ViewT& gVVar,
        const ViewT& gVMin, const ViewT& gVMax,
        // RadialProjectionSingleQuantity device views (binBounds)
        const ViewT& gBAvg, const ViewT& gBVar,
        const ViewT& gBMin, const ViewT& gBMax,
        // Multiplicities (may be a different scalar type)
        const ViewM& gMult,
        ptrdiff_t b) const
    {
      if (mult[b] == T(0)) return;  // No contributions — skip to preserve identity values in global arrays

      Kokkos::atomic_add(&gVAvg(b), vAvg[b]);
      Kokkos::atomic_add(&gVVar(b), vVar[b]);
      Kokkos::atomic_min(&gVMin(b), vMin[b]);
      Kokkos::atomic_max(&gVMax(b), vMax[b]);

      Kokkos::atomic_add(&gBAvg(b), bAvg[b]);
      Kokkos::atomic_add(&gBVar(b), bVar[b]);
      Kokkos::atomic_min(&gBMin(b), bMin[b]);
      Kokkos::atomic_max(&gBMax(b), bMax[b]);

      Kokkos::atomic_add(&gMult(b), mult[b]);
    }
  };

} // namespace TempLat

#endif // DEVICE_KOKKOS
#endif // TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_RADIALPROJECTOR_SCRATCH_H
