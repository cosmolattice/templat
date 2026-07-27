#ifndef TEMPLAT_FFT_EXTERNAL_FFTW_FFTWINTERFACE_H
#define TEMPLAT_FFT_EXTERNAL_FFTW_FFTWINTERFACE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include "TempLat/fft/external/fftw/fftwguard.h"
#include "TempLat/fft/external/fftw/fftwmemorylayout.h"
#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/fft/ffttopology.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"

#include <cstdint>
#include <limits>

namespace TempLat
{
  MakeException(FFTWMPISendrecvCountOverflowException);

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

    virtual device::Idx getMaximumNumberOfDimensionsToDivide(device::Idx nDimensions) { return 1; };

    /** @brief FFTW (including fftw-mpi) uses a slab decomposition: split only the leading dimension.
     *
     *  The grid is pinned to [P, 1, ...] rather than left zeroed for `MPI_Dims_create` to fill.
     *  FFTW slabs the leading dimension across the ranks of its communicator in rank order, so
     *  that shape is not a choice — it is what FFTW will do. Reporting it explicitly keeps this
     *  in agreement with `topology()` below, and lets the group cross-check compare exact dims
     *  instead of merely counting how many dimensions are split. */
    static FFTDecomposition<NDim> decomposition([[maybe_unused]] MPICommReference baseComm,
                                                [[maybe_unused]] device::IdxArray<NDim> nGridPoints)
    {
#ifdef HAVE_MPI
      // FFTW-MPI's slab decomposition does its global transpose with pairwise MPI_Sendrecv calls
      // whose count argument is `int`. The per-rank-pair chunk is approximately
      //   (N0 / P) * (N1 / P) * prod(N_i, i >= 2)
      // elements (in MPI_DOUBLE-equivalent units). Once this exceeds INT_MAX, FFTW-MPI fails at
      // runtime inside the FFT with MPI_ERR_ARG ("invalid count argument"). FFTW-MPI does not
      // chunk-split the transpose, so the only fixes are: more ranks, smaller lattice, or use
      // ParaFaFT (pencil decomposition + 64-bit-safe transposes).
      if constexpr (NDim >= 2) {
        const int nProcesses = baseComm.size();
        if (nProcesses >= 2) {
          const uint64_t kIntMax = static_cast<uint64_t>(std::numeric_limits<int>::max());
          auto perPairCount = [&](int p) {
            uint64_t b0 = (static_cast<uint64_t>(nGridPoints[0]) + p - 1) / p;
            uint64_t b1 = (static_cast<uint64_t>(nGridPoints[1]) + p - 1) / p;
            uint64_t inner = 1;
            for (size_t i = 2; i < NDim; ++i)
              inner *= static_cast<uint64_t>(nGridPoints[i]);
            return b0 * b1 * inner;
          };
          const uint64_t perPair = perPairCount(nProcesses);
          if (perPair > kIntMax) {
            int minP = nProcesses + 1;
            for (; minP <= nProcesses * 1024; ++minP)
              if (perPairCount(minP) <= kIntMax) break;
            throw FFTWMPISendrecvCountOverflowException(
                "FFTW-MPI 1-D slab decomposition with ", nProcesses,
                " ranks would require an MPI_Sendrecv with count ~", perPair,
                " elements during its global transpose, exceeding INT_MAX (", kIntMax,
                "). FFTW-MPI uses an `int` count and does not split the transpose, so this fails "
                "at runtime with MPI_ERR_ARG (\"invalid count argument\"). Either (a) rebuild "
                "with ParaFaFT (-DPARAFAFT=ON), which uses a pencil decomposition that avoids "
                "this limit; or (b) use at least ",
                minP, " MPI ranks; or (c) reduce the lattice size.");
          }
        }
      }
#endif
      FFTDecomposition<NDim> result{};
      for (size_t i = 0; i < NDim; ++i)
        result.dims[i] = 1;
      result.dims[0] = baseComm.size();
      result.nDimsToSplit = baseComm.size() > 1 ? 1 : 0;
      return result;
    }

    /** @brief The MPI topology FFTW-MPI will use.
     *
     *  FFTW-MPI builds no Cartesian communicator of its own: it slabs the leading dimension
     *  across the ranks of whatever communicator it is given, in that communicator's rank order.
     *  So the base communicator *is* the authoritative topology, the grid is [P, 1, ...], and
     *  this rank sits at [rank, 0, ...]. Returning it explicitly (rather than leaving `dims`
     *  zeroed for `MPI_Dims_create`) is what lets the caller build a Cartesian group whose
     *  coordinates provably match the slab FFTW will actually produce.
     *
     *  Runs the same INT_MAX sendrecv-overflow check as `decomposition`.
     */
    static FFTTopology<NDim> topology(MPICommReference baseComm, device::IdxArray<NDim> nGridPoints)
    {
      // Also runs the INT_MAX sendrecv-overflow guard.
      const FFTDecomposition<NDim> shape = decomposition(baseComm, nGridPoints);

      FFTTopology<NDim> result{};
      result.comm = baseComm;
      result.dims = shape.dims;
      result.nDimsToSplit = shape.nDimsToSplit;
      for (size_t i = 0; i < NDim; ++i)
        result.coords[i] = 0;
      result.coords[0] = baseComm.rank(); // slab index == rank on the slabbed leading axis

      result.checkInvariants(baseComm.size());
      return result;
    }

    virtual IntrinsicScales getIntrinsicRescaleToGetUnnormalizedFFT(device::Idx nGridPoints) { return {}; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };

} // namespace TempLat

#endif
