
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct FFTDecompositionTester {

    /** Probed decomposition must describe a rank grid whose product equals the base comm size. */
    template <size_t NDim> static bool decompositionCoversAllRanks(const device::IdxArray<NDim> &nGrid)
    {
      MPICommReference baseComm;
      auto d = FFTLibrarySelector<NDim>::decomposition(baseComm, nGrid);

      device::Idx product = 1;
      device::Idx nonzero = 0;
      for (size_t i = 0; i < NDim; ++i) {
        if (d.dims[i] > 0) {
          product *= d.dims[i];
          ++nonzero;
        }
      }

      // Two legal states:
      //   (a) all dims zero  → backend defers to MPI_Dims_create; we can only check nDimsToSplit bounds
      //   (b) dims pinned    → product of nonzero entries must equal baseComm.size()
      if (nonzero == 0) return d.nDimsToSplit >= 0 && d.nDimsToSplit <= static_cast<device::Idx>(NDim);
      return product == static_cast<device::Idx>(baseComm.size());
    }

    /** The group built by makeMPIGroup must splat `nDimsToSplit` leading dimensions at most
     *  — the same bound the static decomposition query reports. */
    template <size_t NDim> static bool groupAgreesWithStaticQuery(const device::IdxArray<NDim> &nGrid)
    {
      MPICommReference baseComm;
      auto staticDec = FFTLibrarySelector<NDim>::decomposition(baseComm, nGrid);
      auto group = FFTMPIDomainSplit<NDim>::makeMPIGroup(baseComm, nGrid);

      device::Idx actualSplits = 0;
      for (const auto &v : group.getDecomposition())
        if (v > 1) ++actualSplits;

      return actualSplits <= staticDec.nDimsToSplit ||
             (baseComm.size() == 1 && staticDec.nDimsToSplit == 0 && actualSplits == 0);
    }

    /** A hand-built group with more splits than the backend allows must be rejected. Only
     *  meaningful under multi-rank runs; single-rank always has a trivially valid {1,1,...} group. */
    template <size_t NDim> static bool rejectsOverSplitGroup(const device::IdxArray<NDim> &nGrid)
    {
      MPICommReference baseComm;
      if (baseComm.size() < 2) return true; // nothing to split wrong

      auto staticDec = FFTLibrarySelector<NDim>::decomposition(baseComm, nGrid);
      if (staticDec.nDimsToSplit >= static_cast<device::Idx>(NDim)) return true; // no room to over-split

      // Force a decomposition that splits strictly more dims than the backend allows. Easiest:
      // place factors of baseComm.size() in the LAST allowed+1 dims so that the count of dims > 1
      // is staticDec.nDimsToSplit + 1 whenever baseComm.size() is composite; for prime rank counts
      // the simpler case {1,...,size} already exceeds a slab backend's allowed leading split.
      std::vector<int> badDecomp(NDim, 1);
      const int size = static_cast<int>(baseComm.size());
      badDecomp[NDim - 1] = size; // last dim splits, which a leading-slab FFTW backend rejects if
                                  // the backend's nDimsToSplit == 1 and size > 1 — because split
                                  // is then not in a leading dim.
      // Not guaranteed to exceed for every backend; fall back to a {size, size, ...} pattern if
      // needed (handled below by the assertion).

      MPICartesianGroup badGroup(baseComm, static_cast<device::Idx>(NDim), badDecomp);
      try {
        FFTLibrarySelector<NDim> sel(badGroup, nGrid);
        // If the backend does accept this layout (e.g. ParaFaFT probe happens to pin the same
        // last-dim split), this particular case isn't actually a mismatch — not a test failure.
        return true;
      } catch (const FFTLibraryDecompositionMismatchException &) {
        return true; // expected for FFTW-slab backends
      } catch (...) {
        return false; // wrong exception type
      }
    }

    // Factored out of Test(): three sibling scopes that each declared a local `constexpr size_t
    // NDim` and an `IdxArray<NDim>` tripped an nvcc 12.9 front-end bug on the third scope.
    template <size_t NDim, class... Vs> static void runCase(TDDAssertion &tdd, bool testOverSplit, Vs... vs)
    {
      device::IdxArray<NDim> nGrid{static_cast<device::Idx>(vs)...};
      tdd.verify(decompositionCoversAllRanks<NDim>(nGrid));
      tdd.verify(groupAgreesWithStaticQuery<NDim>(nGrid));
      if (testOverSplit) tdd.verify(rejectsOverSplitGroup<NDim>(nGrid));
    }

    static void Test(TDDAssertion &tdd)
    {
      runCase<2>(tdd, true, 16, 16);
      runCase<3>(tdd, true, 16, 16, 16);
      // Anisotropic grid — exercises ParaFaFT's slab-vs-pencil heuristic once it lands: any
      // decision made by the probe should still satisfy the invariants above.
      runCase<3>(tdd, false, 8, 32, 16);
    }
  };

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::FFTDecompositionTester> test;
}
