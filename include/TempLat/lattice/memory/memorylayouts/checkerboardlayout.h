#ifndef TEMPLAT_LATTICE_MEMORY_MEMORYLAYOUTS_CHECKERBOARDLAYOUT_H
#define TEMPLAT_LATTICE_MEMORY_MEMORYLAYOUTS_CHECKERBOARDLAYOUT_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#include "TempLat/lattice/memory/memorylayouts/layoutstruct.h"

namespace TempLat
{
  enum class Parity : device::Idx { Even = 0, Odd = 1 };

  /** @brief Wraps a LayoutStruct to iterate over only even or odd parity lattice sites.
   *
   *  Halves the last memory dimension (NDim-1) and provides a reconstruction formula
   *  that maps half-indices back to full memory indices matching the requested parity.
   *
   *  Usage: pass to device::iteration::foreach or device::iteration::reduce.
   *  The functor receives fully reconstructed full-lattice indices.
   *
   * Unit test: ctest -R test-checkerboardlayout
   **/
  template <size_t NDim> struct CheckerboardLayout {

    CheckerboardLayout(const LayoutStruct<NDim> &layout, Parity parity)
    {
      const auto &sizesInMemory = layout.getSizesInMemory();
      const auto &localStarts = layout.getLocalStarts();
      const device::Idx nGhosts = layout.getNGhosts();

      for (size_t d = 0; d < NDim; ++d) {
        mStarts[d] = nGhosts;
        mStops[d] = sizesInMemory[d];
      }
      // Halve the last memory dimension
      mStops[NDim - 1] = (sizesInMemory[NDim - 1] + 1) / 2;

      // Precompute parity constant:
      //   parityBase = (target_parity + sum_g localStarts[g] + (D-1)*nGhosts) mod 2
      device::Idx sumLocalStarts = 0;
      for (size_t g = 0; g < NDim; ++g)
        sumLocalStarts += localStarts[g];

      mParityBase =
          ((static_cast<device::Idx>(parity) + sumLocalStarts + static_cast<device::Idx>(NDim - 1) * nGhosts) % 2 + 2) %
          2;

      mNGhosts = nGhosts;
      mFullSize = sizesInMemory[NDim - 1];
    }

    DEVICE_INLINE_FUNCTION
    const device::IdxArray<NDim> &getStarts() const { return mStarts; }

    DEVICE_INLINE_FUNCTION
    const device::IdxArray<NDim> &getStops() const { return mStops; }

    /** @brief Reconstruct the full memory index from a half-index. */
    DEVICE_INLINE_FUNCTION
    device::IdxArray<NDim> reconstruct(const device::IdxArray<NDim> &halfIdx) const
    {
      device::IdxArray<NDim> fullIdx = halfIdx;

      // Compute offset from parity of the other dimensions
      device::Idx otherSum = 0;
      for (size_t d = 0; d + 1 < NDim; ++d)
        otherSum += halfIdx[d];

      const device::Idx offset = (mParityBase + otherSum) % 2;

      fullIdx[NDim - 1] = mNGhosts + 2 * (halfIdx[NDim - 1] - mNGhosts) + offset;
      return fullIdx;
    }

    /** @brief Check whether a reconstructed full index is within bounds.
     *  Out-of-bounds occurs at most once per row when the local size along
     *  the split dimension is odd.
     */
    DEVICE_INLINE_FUNCTION
    bool isInBounds(const device::IdxArray<NDim> &fullIdx) const { return fullIdx[NDim - 1] < mNGhosts + mFullSize; }

  private:
    device::IdxArray<NDim> mStarts;
    device::IdxArray<NDim> mStops;
    device::Idx mParityBase;
    device::Idx mNGhosts;
    device::Idx mFullSize;
  };

  // ---- Device-side functor wrappers for foreach/reduce ----

  /** @brief Wraps a foreach functor to apply checkerboard reconstruction. */
  template <size_t NDim, typename Functor> struct CheckerboardForEachWrapper {
    CheckerboardLayout<NDim> mCB;
    Functor mFunctor;

    DEVICE_INLINE_FUNCTION
    void operator()(const device::IdxArray<NDim> &halfIdx) const
    {
      const auto fullIdx = mCB.reconstruct(halfIdx);
      if (mCB.isInBounds(fullIdx)) mFunctor(fullIdx);
    }
  };

  /** @brief Wraps a reduce functor to apply checkerboard reconstruction. */
  template <size_t NDim, typename Functor> struct CheckerboardReduceWrapper {
    CheckerboardLayout<NDim> mCB;
    Functor mFunctor;

    template <typename T> DEVICE_INLINE_FUNCTION void operator()(const device::IdxArray<NDim> &halfIdx, T &update) const
    {
      const auto fullIdx = mCB.reconstruct(halfIdx);
      if (mCB.isInBounds(fullIdx)) mFunctor(fullIdx, update);
    }
  };

} // namespace TempLat

#endif
