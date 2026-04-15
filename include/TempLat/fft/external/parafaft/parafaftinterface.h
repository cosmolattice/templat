#if !defined(TempLat_FFT_EXTERNAL_PARAFAFT_PARAFAFTINTERFACE_H) && defined(HAVE_PARAFAFT)
#define TempLat_FFT_EXTERNAL_PARAFAFT_PARAFAFTINTERFACE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#include "TempLat/fft/external/parafaft/parafaftmemorylayout.h"
#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"

namespace TempLat
{
  /** @brief Top-level interface for parafaft FFT backend.
   *
   * Parafaft uses pencil decomposition to parallelize FFTs across D-1 dimensions.
   * For 3D, this means 2D parallelization, but using MPI_Alltoallw
   * instead of local transposes.
   *
   * Inheritance hierarchy:
   * ParafaftInterface -> ParafaftMemoryLayout -> ParafaftPlanner -> FFTLibraryInterface
   *
   * Unit test: ctest -R test-parafaftinterface
   **/
  template <size_t NDim> class ParafaftInterface : public ParafaftMemoryLayout<NDim>
  {
  public:
    ParafaftInterface() {}

    /**
     * @brief Returns maximum number of dimensions that can be parallelized.
     *
     * Parafaft supports D-1 dimensional decomposition (pencil decomposition).
     * For 3D data, this returns 2 (distribute over 2 dimensions).
     */
    virtual ptrdiff_t getMaximumNumberOfDimensionsToDivide(ptrdiff_t nDimensions) override
    {
      // Pencil decomposition: all but one dimension can be distributed
      return std::max((ptrdiff_t)1, nDimensions - 1);
    }

    /** @brief Pencil decomposition for ParaFaFT, probed from a temporary `ParaFaFT_R2C` planner.
     *
     *  ParaFaFT's own heuristic (which may pick slab or pencil depending on rank count and
     *  lattice extents in the current/in-progress extension) is the source of truth. We build
     *  a throwaway planner on the base comm and query `get_domain_decomposition` to lock in
     *  the exact grid shape that ParaFaFT will use later on the real comm. Because ParaFaFT's
     *  choice is a deterministic function of `(baseComm.size(), nGridPoints)`, the same probe
     *  on every rank yields identical dims.
     */
    static FFTDecomposition<NDim> decomposition(MPICommReference baseComm, device::IdxArray<NDim> nGridPoints)
    {
      FFTDecomposition<NDim> result{};
      int globalShape[NDim];
      for (size_t i = 0; i < NDim; ++i) globalShape[i] = static_cast<int>(nGridPoints[i]);

      parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend> probe(globalShape, baseComm);

      int probeDims[NDim];
      probe.get_domain_decomposition(probeDims);

      ptrdiff_t nSplit = 0;
      for (size_t i = 0; i < NDim; ++i) {
        result.dims[i] = probeDims[i];
        if (probeDims[i] > 1) ++nSplit;
      }
      result.nDimsToSplit = nSplit;
      return result;
    }

    /**
     * @brief Returns intrinsic rescaling factors.
     *
     * Parafaft (via FFTW) produces unnormalized FFTs, same as FFTW.
     * Returns default IntrinsicScales (1.0, 1.0).
     */
    virtual IntrinsicScales getIntrinsicRescaleToGetUnnormalizedFFT(ptrdiff_t nGridPoints) override
    {
      return IntrinsicScales();
    }
  };
} // namespace TempLat

#endif
