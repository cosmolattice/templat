
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/util/tdd/tdd.h"

#include <numeric>
#include <vector>

namespace TempLat
{
  /** @brief Tests that the per-rank layouts ParaFaFT reports actually tile the global lattice.
   *
   *  Every rank believing it owns a plausible-looking subdomain is not enough: the subdomains
   *  must jointly cover the lattice exactly once. An ordering disagreement between TempLat's
   *  communicator and ParaFaFT's produces per-rank layouts that each look fine while overlapping
   *  and leaving holes globally — which is why this checks the union, not the parts.
   *
   *  Only meaningful when ParaFaFT is the selected backend; degrades to a pass otherwise.
   */
  struct ParafaftMemoryLayoutTester {

#if defined(HAVE_MPI) && defined(HAVE_PARAFAFT)
    /** Gather every rank's [start, start+size) box and check it tiles the global volume once. */
    template <size_t NDim>
    static bool tilesCoverLatticeExactlyOnce(const device::IdxArray<NDim> &nGrid, bool fourierSpace)
    {
      MPICommReference baseComm;
      auto group = FFTMPIDomainSplit<NDim>::makeMPIGroup(baseComm, nGrid);
      FFTLibrarySelector<NDim> selector(group, nGrid);
      const auto &layout = selector.getLayout();
      const auto &space = fourierSpace ? layout.fourierSpace : layout.configurationSpace;

      // Global extents: the Fourier last axis is the halved r2c one.
      std::vector<long> global(NDim);
      for (size_t i = 0; i < NDim; ++i)
        global[i] = static_cast<long>(space.getGlobalSizes()[i]);

      const int nRanks = static_cast<int>(baseComm.size());
      std::vector<int> myBox(2 * NDim);
      for (size_t i = 0; i < NDim; ++i) {
        myBox[i] = static_cast<int>(space.getLocalStarts()[i]);
        myBox[NDim + i] = static_cast<int>(space.getLocalSizes()[i]);
      }

      std::vector<int> allBoxes(2 * NDim * nRanks);
      MPI_Allgather(myBox.data(), 2 * NDim, MPI_INT, allBoxes.data(), 2 * NDim, MPI_INT, baseComm);

      // Count how many ranks claim each site. Every site must be claimed exactly once.
      long totalSites = std::accumulate(global.begin(), global.end(), 1L, std::multiplies<long>());
      if (totalSites <= 0 || totalSites > (1L << 24)) return false; // keep the test cheap
      std::vector<int> claims(totalSites, 0);

      for (int r = 0; r < nRanks; ++r) {
        const int *start = &allBoxes[2 * NDim * r];
        const int *size = start + NDim;

        long boxSites = 1;
        for (size_t i = 0; i < NDim; ++i)
          boxSites *= size[i];

        // Walk the box in row-major order and mark each site it claims.
        std::vector<long> coord(NDim);
        for (long flat = 0; flat < boxSites; ++flat) {
          long rem = flat;
          for (size_t i = NDim; i-- > 0;) {
            coord[i] = start[i] + (rem % size[i]);
            rem /= size[i];
            if (coord[i] < 0 || coord[i] >= global[i]) return false; // box escapes the lattice
          }

          long globalFlat = 0;
          for (size_t i = 0; i < NDim; ++i)
            globalFlat = globalFlat * global[i] + coord[i];

          ++claims[globalFlat];
        }
      }

      for (long s = 0; s < totalSites; ++s)
        if (claims[s] != 1) return false;

      return true;
    }

    /** The last axis must never be distributed; three unguarded call sites depend on it. */
    template <size_t NDim> static bool lastAxisStaysLocal(const device::IdxArray<NDim> &nGrid)
    {
      MPICommReference baseComm;
      const auto topology = FFTLibrarySelector<NDim>::topology(baseComm, nGrid);
      return topology.dims[NDim - 1] == 1 && topology.coords[NDim - 1] == 0;
    }
#endif

    template <size_t NDim, class... Vs> static void runCase([[maybe_unused]] TDDAssertion &tdd, Vs... vs)
    {
#if defined(HAVE_MPI) && defined(HAVE_PARAFAFT)
      device::IdxArray<NDim> nGrid{static_cast<device::Idx>(vs)...};
      tdd.verify(lastAxisStaysLocal<NDim>(nGrid));
      tdd.verify(tilesCoverLatticeExactlyOnce<NDim>(nGrid, /*fourierSpace=*/false));
      tdd.verify(tilesCoverLatticeExactlyOnce<NDim>(nGrid, /*fourierSpace=*/true));
#else
      ((void)vs, ...);
#endif
    }

    static void Test(TDDAssertion &tdd)
    {
      runCase<2>(tdd, 16, 16);
      runCase<3>(tdd, 16, 16, 16);
      runCase<3>(tdd, 8, 32, 16); // anisotropic: exercises the slab-vs-pencil heuristic
    }
  };

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ParafaftMemoryLayoutTester> test;
}
