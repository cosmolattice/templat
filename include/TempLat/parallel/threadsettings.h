#ifndef TEMPLAT_PARALLEL_THREADS_THREADSETTINGS_H
#define TEMPLAT_PARALLEL_THREADS_THREADSETTINGS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler,  Year: 2025

#include <thread>
#include <algorithm>

#include "TempLat/parallel/device.h"

#include "TempLat/util/log/saycomplete.h"

namespace TempLat
{
  /** @brief A class which does bookkeeping, number of threads allowed by hardware, number of threads allowed by user,
   *etc.
   *
   * Unit test: ctest -R test-threadsettings
   **/
  class ThreadSettings
  {
  public:
    // Put public methods here. These should change very little over time.

    /** @brief Yes, a global variable. Inline, so it remains header only. Constructor is private.
        All you can do with this instance is pass it to a stream:
        `std::cout << ThreadSettings::getInstance() << "\n";`
    */
    static inline ThreadSettings &getInstance()
    {
      static ThreadSettings store;
      return store;
    }

    static void setMPILocalSize(device::Idx newSize) { getInstance().pSetMPILocalSize(newSize); }

    static device::Idx getMPILocalSize() { return getInstance().pGetMPILocalSize(); }

    static void setMPIThreadsNotOK() { getInstance().pSetMPIThreadsNotOK(); }

    static device::Idx getMaxThreadCount() { return getInstance().pGetMaxThreadCount(); }

    friend std::ostream &operator<<(std::ostream &stream, const ThreadSettings &fts)
    {
      stream << "Threading information:\n";
      stream << " - number of cores on machine:                               " << fts.mHardwareNumCores << "\n";
      stream << " - number of mpi processes for this session on this machine: " << fts.mMPILocalSize << "\n";
      stream << " - user specified number of threads per process:             " << fts.mUserAllowedThreadsPerProcess
             << "\n";
      stream << " - resulting maximum number of threads per process:          " << fts.mHardwareAllowedThreadsPerProcess
             << "\n";
      return stream;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    device::Idx mMPILocalSize;
    device::Idx mHardwareNumCores;
    device::Idx mHardwareAllowedThreadsPerProcess;
    device::Idx mUserAllowedThreadsPerProcess;

    ThreadSettings()
        : mMPILocalSize(1), mHardwareNumCores(std::thread::hardware_concurrency()),
          mHardwareAllowedThreadsPerProcess(mHardwareNumCores / mMPILocalSize),
          mUserAllowedThreadsPerProcess(mHardwareAllowedThreadsPerProcess)
    {
      mHardwareAllowedThreadsPerProcess = std::max(mHardwareAllowedThreadsPerProcess, device::Idx{1});

      device::Idx kokkosThreads = 0;
      if (const char *env_p = std::getenv("KOKKOS_NUM_THREADS")) kokkosThreads = std::stoi(env_p);

      device::Idx ompThreads = 0;
      if (const char *env_p = std::getenv("OMP_NUM_THREADS")) ompThreads = std::stoi(env_p);

      mUserAllowedThreadsPerProcess = std::max(kokkosThreads, ompThreads);
      if (kokkosThreads > 0 && ompThreads > 0 && kokkosThreads != ompThreads) {
        sayShort << "Warning: both KOKKOS_NUM_THREADS and OMP_NUM_THREADS are set, but to different values. Using the "
                    "largest value.\n    KOKKOS_NUM_THREADS = "
                 << kokkosThreads << ", OMP_NUM_THREADS = " << ompThreads << ".\n";
      }

      if (mUserAllowedThreadsPerProcess > 0)
        mHardwareAllowedThreadsPerProcess = std::min(mUserAllowedThreadsPerProcess, mHardwareAllowedThreadsPerProcess);

      // std::cout << "Result: mHardwareNumCores = " << mHardwareNumCores << ", mMPILocalSize = " << mMPILocalSize
      //           << ", mUserAllowedThreadsPerProcess = " << mUserAllowedThreadsPerProcess
      //           << ", mHardwareAllowedThreadsPerProcess = " << mHardwareAllowedThreadsPerProcess << "\n";
    }

    void pSetMPILocalSize(device::Idx newSize)
    {
      mMPILocalSize = newSize > 0 ? newSize : 1;
      mHardwareAllowedThreadsPerProcess = std::max(mHardwareNumCores / mMPILocalSize, device::Idx{1});

      if (mUserAllowedThreadsPerProcess > 0)
        mHardwareAllowedThreadsPerProcess = std::min(mUserAllowedThreadsPerProcess, mHardwareAllowedThreadsPerProcess);
    }

    device::Idx pGetMPILocalSize() const { return mMPILocalSize; }

    void pSetMPIThreadsNotOK() { mUserAllowedThreadsPerProcess = 1; }

    device::Idx pGetMaxThreadCount() const { return mHardwareAllowedThreadsPerProcess; }
  };
} // namespace TempLat

#endif
