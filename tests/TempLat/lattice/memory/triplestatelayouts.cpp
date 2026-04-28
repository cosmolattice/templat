
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019
#include "TempLat/lattice/memory/triplestatelayouts.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  template <size_t NDim> struct TripleStateLayoutsTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void TripleStateLayoutsTester<NDim>::Test(TDDAssertion &tdd)
  {
    constexpr device::Idx nGhost = 10;
    device::IdxArray<4> nGrid{256, 256, 256, 256};
    auto group = FFTMPIDomainSplit<4>::makeMPIGroup(MPICommReference(), nGrid);
    FFTLibrarySelector<4> theLibrary(group, nGrid);

    TripleStateLayouts<4> tsl(theLibrary.getLayout(), nGhost);

    if (TDDRegister::isSingleUnitTest()) {
      say << tsl << "\n";
    }

    /* Compute expected memory from the per-rank local sizes the backend chose,
     * so this works for any decomposition (slab, pencil, ...). */
    device::Idx expectedGhostMem = 1;
    for (auto&& s : theLibrary.getLayout().configurationSpace.getLocalSizes())
      expectedGhostMem *= (s + 2 * nGhost);
    device::Idx expectedFFTMem = theLibrary.getLayout().getMinimalMemorySize();
    device::Idx expectedMem = std::max(expectedGhostMem, expectedFFTMem);

    tdd.verify(tsl.getNecessaryMemoryAllocation() == expectedMem);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::TripleStateLayoutsTester<3>> test;
}
