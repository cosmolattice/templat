#if !defined(TempLat_FFT_EXTERNAL_KOKKOSFFT_KOKKOSFFTINTERFACE_H) && defined(HAVE_KOKKOSFFT)
#define TEMPLAT_FFT_EXTERNAL_KOKKOSFFT_KOKKOSFFTINTERFACE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2025

#include "TempLat/fft/external/kokkosfft/kokkosfftmemorylayout.h"
#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/fft/ffttopology.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"

namespace TempLat
{
  /** @brief A class which implements all of FFTLibraryInterface. The larger methods are implemented in classes from
   *which we inherit, in a linear chain: KokkosFFTMemoryLayout and KokkosFFTPlanner.
   *
   *
   * Unit test: ctest -R test-kokkosinterface
   **/
  template <size_t NDim> class KokkosFFTInterface : public KokkosFFTMemoryLayout<NDim>
  {
  public:
    // Put public methods here. These should change very little over time.
    KokkosFFTInterface() = default;

    /** @brief KokkosFFT does not run under MPI, so it can divide nothing.
     *
     *  This used to return `nDimensions - 1`, contradicting its own `decomposition()` report of
     *  zero split dimensions.
     */
    virtual device::Idx getMaximumNumberOfDimensionsToDivide(device::Idx /*nDimensions*/) { return 0; };

    /** @brief KokkosFFT does not run under MPI: the entire lattice is local. Report "no split".
     *  The selector only reaches this backend when `baseComm.size() == 1`. */
    static FFTDecomposition<NDim> decomposition(MPICommReference /*baseComm*/, device::IdxArray<NDim> /*nGridPoints*/)
    {
      return FFTDecomposition<NDim>{/*nDimsToSplit=*/0, {}};
    }

    /** @brief Trivial topology: one cell, this rank at the origin, on the base communicator.
     *  The selector only reaches this backend at a single rank, so the grid is all ones. */
    static FFTTopology<NDim> topology(MPICommReference baseComm, device::IdxArray<NDim> /*nGridPoints*/)
    {
      FFTTopology<NDim> result{};
      result.comm = baseComm;
      for (size_t i = 0; i < NDim; ++i) {
        result.dims[i] = 1;
        result.coords[i] = 0;
      }
      result.nDimsToSplit = 0;

      result.checkInvariants(baseComm.size());
      return result;
    }

    virtual IntrinsicScales getIntrinsicRescaleToGetUnnormalizedFFT(device::Idx nGridPoints)
    {
      return IntrinsicScales{1., 1};
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };

} // namespace TempLat

#endif
