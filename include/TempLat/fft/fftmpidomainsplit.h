#ifndef TEMPLAT_FFT_FFTMPIDOMAINSPLIT_H
#define TEMPLAT_FFT_FFTMPIDOMAINSPLIT_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019

#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesiangroup.h"

namespace TempLat
{
  /** @brief A class which combines the MPIDomainSplit with the FFTLibrarySelector limits on domain splitting.
   *
   * Unit test: ctest -R test-fftmpidomainsplit
   **/

  template <size_t NDim> class FFTMPIDomainSplit
  {
  public:
    // Put public methods here. These should change very little over time.

    /** @brief Build the MPI domain decomposition that matches the FFT backend selected for this
     *  `(baseGroup.size(), nGridPoints)` combination. If the backend pinned an explicit grid shape
     *  (e.g. ParaFaFT's probe), that shape is honoured verbatim; otherwise `MPIDomainSplit` factors
     *  the rank count over the backend's allowed number of split dimensions.
     */
    static std::vector<int> makeDomainDecomposition(MPICommReference baseGroup,
                                                    const device::IdxArray<NDim> &nGridPoints)
    {
      FFTDecomposition<NDim> d = FFTLibrarySelector<NDim>::decomposition(baseGroup, nGridPoints);

      bool explicitGrid = d.nDimsToSplit > 0;
      for (device::Idx i = 0; i < d.nDimsToSplit; ++i)
        if (d.dims[i] <= 0) {
          explicitGrid = false;
          break;
        }

      if (explicitGrid) {
        std::vector<int> out(NDim, 1);
        for (size_t i = 0; i < NDim; ++i)
          out[i] = std::max(1, d.dims[i]);
        return out;
      }

      MPIDomainSplit theSplit(baseGroup.size(), static_cast<device::Idx>(NDim), d.nDimsToSplit);
      return theSplit;
    }

    static MPICartesianGroup makeMPIGroup(MPICommReference baseGroup, const device::IdxArray<NDim> &nGridPoints)
    {
      return MPICartesianGroup(baseGroup, static_cast<device::Idx>(NDim),
                               makeDomainDecomposition(baseGroup, nGridPoints));
    }
    /* default using comm_world */
    static MPICartesianGroup makeMPIGroup(const device::IdxArray<NDim> &nGridPoints)
    {
      return makeMPIGroup(MPICommReference(), nGridPoints);
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    FFTMPIDomainSplit() = default;
  };

} // namespace TempLat

#endif
