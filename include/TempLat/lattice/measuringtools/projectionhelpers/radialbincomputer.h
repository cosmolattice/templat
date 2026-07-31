#ifndef TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_RADIALBINCOMPUTER_H
#define TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_RADIALBINCOMPUTER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include <algorithm>

#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A class which returns an integer bin for a given fp value.
   *
   * Unit test: ctest -R test-radialbincomputer
   **/
  template <typename T> class RadialBinComputer
  {
  public:
    // Templated on the projected quantity's scalar type: with hard-coded doubles, the binning kernel
    // performed an fp64 divide and floor per lattice site even for a single-precision field.
    RadialBinComputer(T minVal, T maxVal, device::Idx nBins, T deltaKBin)
        : mMinVal(minVal), mMaxVal(maxVal), mRange(mMaxVal - mMinVal), mNBins(nBins), mHighestBin(nBins - 1),
          mDeltakBin(deltaKBin)
    {
      if (mRange <= 0) mRange = 1;
    }

    /** @brief Call this for your value, receive a bin index in return. */
    DEVICE_FUNCTION
    device::Idx operator()(T value) const
    {
      const device::Idx bin = static_cast<device::Idx>(device::floor((value - mMinVal) / mDeltakBin));
      return device::min(mHighestBin, device::max(device::Idx(0), bin));
    }

    template <typename S> void setCentralBinBounds(std::vector<S> &res)
    {
      res = std::vector<S>(mNBins);
      const T steps = mDeltakBin;
      for (device::Idx i = 0; i < mNBins; ++i) {
        res[i] = mMinVal + mDeltakBin / 2 + i * steps;
      }
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    T mMinVal;
    T mMaxVal;
    T mRange;
    device::Idx mNBins;
    device::Idx mHighestBin;
    T mDeltakBin;
  };
} // namespace TempLat

#endif
