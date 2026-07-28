
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/memory/memorymanager.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  template <typename T, size_t NDim> struct MemoryManagerTester {
    static void Test(TDDAssertion &tdd);
  };

  template <typename T, size_t NDim> void MemoryManagerTester<T, NDim>::Test(TDDAssertion &tdd)
  {
    auto toolBox = MemoryToolBox<NDim>::makeShared(192, 2);

    toolBox->mFFTLibrary.setVerbose();

    MemoryManager<T, NDim> mManager(toolBox);

    if (TDDRegister::isSingleUnitTest()) {
      std::cerr << mManager << "\n\n";

      toolBox->setVerbose();
    }

    /* first allocation */
    tdd.verify(mManager.confirmConfigSpace() > 0);

    /* fft necessary */
    tdd.verify(mManager.confirmFourierSpace() > 0);

    /* fft not necessary */
    tdd.verify(mManager.confirmFourierSpace() == 0);

    /* fft necessary */
    tdd.verify(mManager.confirmConfigSpace() > 0);

    /* fft not necessary */
    tdd.verify(mManager.confirmConfigSpace() == 0);

    /* ghost update necessary */
    tdd.verify(mManager.confirmGhostsUpToDate() > 0);

    /* ghost update not necessary */
    tdd.verify(mManager.confirmGhostsUpToDate() == 0);

    /* fft not necessary */
    tdd.verify(mManager.confirmConfigSpace() == 0);

    /* ghost update not necessary */
    tdd.verify(mManager.confirmGhostsUpToDate() == 0);

    /* fft necessary */
    tdd.verify(mManager.confirmFourierSpace() > 0);

    /* fft && ghost update necessary */
    tdd.verify(mManager.confirmGhostsUpToDate() > 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MemoryManagerTester<double, 3>> test;
}
