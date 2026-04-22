#ifndef TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_RADIALBINCOMPUTER_H
#define TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_RADIALBINCOMPUTER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019

#include <algorithm>

#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A class which returns an integer bin for a given fp value.
   *
   * Unit test: ctest -R test-radialbincomputer
   **/
  class RadialBinComputer
  {
  public:
    RadialBinComputer(double minVal, double maxVal, device::Idx nBins, double deltaKBin)
        : mMinVal(minVal), mMaxVal(maxVal), mRange(mMaxVal - mMinVal), mNBins(nBins), mHighestBin(nBins - 1),
          mDeltakBin(deltaKBin)
    {
      if (mRange <= 0) mRange = 1;
    }

    /** @brief Call this for your value, receive a bin index in return. */
    DEVICE_FUNCTION
    device::Idx operator()(double value) const
    {
      const device::Idx bin = static_cast<device::Idx>(device::floor((value - mMinVal) / mDeltakBin));
      return device::min(mHighestBin, device::max(device::Idx(0), bin));
    }

    template <typename T> void setCentralBinBounds(std::vector<T> &res)
    {
      res = std::vector<T>(mNBins);
      T steps = mDeltakBin;
      for (device::Idx i = 0; i < mNBins; ++i) {
        res[i] = mMinVal + mDeltakBin / 2. + i * steps;
      }
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    double mMinVal;
    double mMaxVal;
    double mRange;
    device::Idx mNBins;
    device::Idx mHighestBin;
    double mDeltakBin;
  };
} // namespace TempLat

#endif
