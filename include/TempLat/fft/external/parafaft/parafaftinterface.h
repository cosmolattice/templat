#if !defined(TempLat_FFT_EXTERNAL_PARAFAFT_PARAFAFTINTERFACE_H) && defined(HAVE_PARAFAFT)
#define TempLat_FFT_EXTERNAL_PARAFAFT_PARAFAFTINTERFACE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026

#include "TempLat/fft/external/parafaft/parafaftmemorylayout.h"
#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/fft/ffttopology.h"
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
    virtual device::Idx getMaximumNumberOfDimensionsToDivide(device::Idx nDimensions) override
    {
      // Pencil decomposition: all but one dimension can be distributed
      return std::max((device::Idx)1, nDimensions - 1);
    }

    /** @brief The full MPI topology ParaFaFT will use, probed from a temporary `ParaFaFT_R2C`.
     *
     *  ParaFaFT's own heuristic is the source of truth, and we take not just its grid *shape* but
     *  the communicator whose rank order defines it. That matters because ParaFaFT derives
     *  `real_global_start_` / `complex_global_start_` from its own Cartesian coordinates, and
     *  those become TempLat's local starts. Taking the shape alone (as this used to) left
     *  TempLat free to build a differently-ordered communicator whose ghost neighbours did not
     *  match those starts.
     *
     *  The returned communicator is a duplicate that we own; `MPICommReference` refcounts it and
     *  frees it when the last reference dies. It must be a duplicate because ParaFaFT's
     *  destructor frees its internal Cartesian communicator and `probe` is a temporary.
     */
    static FFTTopology<NDim> topology(MPICommReference baseComm, device::IdxArray<NDim> nGridPoints)
    {
      FFTTopology<NDim> result{};
      int globalShape[NDim];
      for (size_t i = 0; i < NDim; ++i)
        globalShape[i] = static_cast<int>(nGridPoints[i]);

      // Pin the probe to double: the decomposition depends only on (baseComm.size(),
      // nGridPoints), not on the floating-point precision, and this static method is called
      // before any T is known. Double is always available; float depends on
      // PARAFAFT_FFTW3F_AVAILABLE.
      parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend<double>> probe(globalShape, baseComm);

      int probeDims[NDim];
      probe.get_domain_decomposition(probeDims);

      // ParaFaFT's grid spans NDim-1 dimensions; the last axis is never distributed, so it has
      // no coordinate of its own and stays 0. Ask ParaFaFT for the grid rank rather than
      // assuming the offset.
      constexpr int gridNDims = parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend<double>>::get_grid_ndims();
      int probeCoords[gridNDims];
      probe.get_grid_coords(probeCoords);

      int nSplit = 0;
      for (size_t i = 0; i < NDim; ++i) {
        result.dims[i] = probeDims[i];
        result.coords[i] = (static_cast<int>(i) < gridNDims) ? probeCoords[i] : 0;
        if (probeDims[i] > 1) ++nSplit;
      }
      result.nDimsToSplit = nSplit;
      result.comm = MPICommReference(probe.dup_cartesian_comm());

      result.checkInvariants(baseComm.size());
      return result;
    }

    /** @brief Shape-only query, for callers that do not need the communicator.
     *
     *  Deliberately does not go through `topology()`: that duplicates a communicator, and this is
     *  called on paths (such as verifying a hand-built group) that only ever look at the shape.
     */
    static FFTDecomposition<NDim> decomposition(MPICommReference baseComm, device::IdxArray<NDim> nGridPoints)
    {
      FFTDecomposition<NDim> result{};
      int globalShape[NDim];
      for (size_t i = 0; i < NDim; ++i)
        globalShape[i] = static_cast<int>(nGridPoints[i]);

      parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend<double>> probe(globalShape, baseComm);

      int probeDims[NDim];
      probe.get_domain_decomposition(probeDims);

      int nSplit = 0;
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
    virtual IntrinsicScales getIntrinsicRescaleToGetUnnormalizedFFT(device::Idx nGridPoints) override
    {
      return IntrinsicScales();
    }
  };
} // namespace TempLat

#endif
