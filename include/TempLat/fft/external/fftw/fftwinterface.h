#ifndef TEMPLAT_FFT_EXTERNAL_FFTW_FFTWINTERFACE_H
#define TEMPLAT_FFT_EXTERNAL_FFTW_FFTWINTERFACE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019

#include "TempLat/fft/external/fftw/fftwguard.h"
#include "TempLat/fft/external/fftw/fftwmemorylayout.h"
#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"

namespace TempLat
{
  /** @brief A class which implements all of FFTLibraryInterface. The larger methods are implemented in classes from
   *which we inherit, in a linear chain: FFTWMemoryLayout and FFTWPlanner.
   *
   *
   * Unit test: ctest -R test-fftwinterface
   **/
  template <size_t NDim> class FFTWInterface : public FFTWMemoryLayout<NDim>
  {
  public:
    // Put public methods here. These should change very little over time.
    FFTWInterface()
    {
      // Sanity check: our complex type should be perfectly compatible with FFTW complex.
      static_assert(sizeof(fftwf_complex) == sizeof(complex<float>));
      static_assert(alignof(fftwf_complex) <= alignof(complex<float>));
      static_assert(sizeof(fftw_complex) == sizeof(complex<double>));
      static_assert(alignof(fftw_complex) <= alignof(complex<double>));
    }

    virtual ptrdiff_t getMaximumNumberOfDimensionsToDivide(ptrdiff_t nDimensions) { return 1; };

    /** @brief FFTW (including fftw-mpi) uses a slab decomposition: split only the leading dimension.
     *  `dims` is left as zeros so MPI_Dims_create can fill the single nonzero slot. */
    static FFTDecomposition<NDim> decomposition(MPICommReference /*baseComm*/,
                                                device::IdxArray<NDim> /*nGridPoints*/)
    {
      return FFTDecomposition<NDim>{/*nDimsToSplit=*/1, {}};
    }

    virtual IntrinsicScales getIntrinsicRescaleToGetUnnormalizedFFT(ptrdiff_t nGridPoints) { return {}; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };

} // namespace TempLat

#endif
