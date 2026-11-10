#ifndef TEMPLAT_LATTICE_MEMORY_SWEEPBUFFER_H
#define TEMPLAT_LATTICE_MEMORY_SWEEPBUFFER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <cstddef>

#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"
#include "TempLat/util/log/saycomplete.h"

namespace TempLat
{
  /** @brief Device scratch memory for assignments that read their own target at neighbouring sites, see
   *  sweepAssign. One per MemoryToolBox, allocated on first use and grown when a larger one is needed. The first
   *  allocation prints a warning on rank 0.
   *
   *  Unit test: ctest -R test-sweepassignment
   **/
  class SweepBuffer
  {
  public:
    /** @brief At least count elements of T in device memory. Valid until the next call. */
    template <typename T> T *get(size_t count)
    {
      const size_t bytes = count * sizeof(T);
      if (bytes > mData.size()) {
        if (!isAllocated() && MPICommReference::worldRank() == 0)
          sayShort << "TempLat warning: an assignment reads its own target at neighbouring sites (e.g. "
                      "phi = LatLapl(phi)). It is evaluated plane by plane, using a scratch buffer of "
                   << (bytes + 1023) / 1024 << " kB. Printed once per MemoryToolBox.\n";
        mData = device::memory::NDView<char, 1>("SweepBuffer", bytes);
      }
      return reinterpret_cast<T *>(mData.data());
    }

    bool isAllocated() const { return mData.size() > 0; }

    /** @brief Number of planes evaluated per sweep step; 0 (the default) chooses it from the device. For tests. */
    void setBandWidth(device::Idx bandWidth) { mBandWidth = bandWidth; }
    device::Idx getBandWidth() const { return mBandWidth; }

  private:
    device::memory::NDView<char, 1> mData;
    device::Idx mBandWidth = 0;
  };
} // namespace TempLat

#endif
