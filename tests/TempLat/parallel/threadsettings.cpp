/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/parallel/threadsettings.h"

namespace TempLat
{
  struct ThreadSettingsTester {
    static void Test(TDDAssertion &tdd);
  };

  void ThreadSettingsTester::Test(TDDAssertion &tdd)
  {
    device::Idx initialThreadCount = ThreadSettings::getMaxThreadCount();

    device::Idx initialMPISize = ThreadSettings::getMPILocalSize();

    ThreadSettings::setMPILocalSize(initialMPISize * 2);

    device::Idx newThreadCount = ThreadSettings::getMaxThreadCount();

    device::Idx manuallyComputedThreadCount = std::max((device::Idx)1, initialThreadCount / 2);

    tdd.verify(manuallyComputedThreadCount <= newThreadCount);

    char *env_p = std::getenv("OMP_NUM_THREADS");
    if (env_p != nullptr) {
      device::Idx ompThreads = std::stoi(env_p);
      tdd.verify(ompThreads >= newThreadCount);
    }

    char *kokkosEnv_p = std::getenv("KOKKOS_NUM_THREADS");
    if (kokkosEnv_p != nullptr) {
      device::Idx kokkosThreads = std::stoi(kokkosEnv_p);
      tdd.verify(kokkosThreads >= newThreadCount);
    }

    ThreadSettings::setMPILocalSize(initialMPISize);

    say << ThreadSettings::getInstance();
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ThreadSettingsTester> test;
}
