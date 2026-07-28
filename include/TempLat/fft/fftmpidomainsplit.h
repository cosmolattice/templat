#ifndef TEMPLAT_FFT_FFTMPIDOMAINSPLIT_H
#define TEMPLAT_FFT_FFTMPIDOMAINSPLIT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

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

    /** @brief Build the Cartesian group on the FFT backend's own communicator.
     *
     *  The group is created over `topology.comm` with no reordering, and its coordinates are
     *  cross-checked against the backend's. So the group does not merely have the same *shape*
     *  as the backend's decomposition — it has the backend's actual rank ordering, which is what
     *  the layout's local starts are derived from.
     */
    static MPICartesianGroup makeMPIGroup(MPICommReference baseGroup, const device::IdxArray<NDim> &nGridPoints)
    {
      const FFTTopology<NDim> topology = FFTLibrarySelector<NDim>::topology(baseGroup, nGridPoints);
      return MPICartesianGroup(topology.comm, static_cast<device::Idx>(NDim), topology.dimsVector(),
                               std::vector<int>(topology.coords.begin(), topology.coords.end()));
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
