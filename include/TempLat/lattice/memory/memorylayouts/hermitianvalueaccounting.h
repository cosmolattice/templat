#ifndef TEMPLAT_LATTICE_MEMORY_MEMORYLAYOUTS_HERMITIANVALUEACCOUNTING_H
#define TEMPLAT_LATTICE_MEMORY_MEMORYLAYOUTS_HERMITIANVALUEACCOUNTING_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019

#include <cstddef>
#include <ostream>
#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A number of unique / independent real and imaginary floating point values in a memory layout.
   *
   *
   * Unit test: ctest -R test-hermitianvalueaccounting
   **/
  class HermitianValueAccounting
  {
  public:
    // Put public methods here. These should change very little over time.
    DEVICE_FUNCTION
    HermitianValueAccounting(device::Idx realCount, device::Idx imaginaryCount)
        : mRealCount(realCount), mImaginaryCount(imaginaryCount)
    {
    }

    DEVICE_INLINE_FUNCTION
    device::Idx getRealValueCount() const { return mRealCount; }

    DEVICE_INLINE_FUNCTION
    device::Idx getImaginaryValueCount() const { return mImaginaryCount; }

    friend std::ostream &operator<<(std::ostream &ostream, const HermitianValueAccounting &hva)
    {
      ostream << hva.mRealCount << " independent real values,\n" << hva.mImaginaryCount << " imaginary values.";
      return ostream;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    device::Idx mRealCount;
    device::Idx mImaginaryCount;
  };

} // namespace TempLat

#endif
