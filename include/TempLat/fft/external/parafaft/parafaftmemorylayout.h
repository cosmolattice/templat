#ifndef TEMPLAT_FFT_EXTERNAL_PARAFAFT_PARAFAFTMEMORYLAYOUT_H
#define TEMPLAT_FFT_EXTERNAL_PARAFAFT_PARAFAFTMEMORYLAYOUT_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#ifndef NOFFT
#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
#include <parafaft_r2c.hpp>
#endif
#endif
#endif

#include "TempLat/fft/external/parafaft/parafaftplanner.h"
#include "TempLat/fft/external/fftw/fftwhermitianpartners.h"
#include "TempLat/lattice/memory/memorylayouts/fftlayoutstruct.h"

#include <numeric>

namespace TempLat
{

  MakeException(ParafaftMemoryLayoutException);

  /** @brief Memory layout computation for parafaft.
   *
   * Parafaft uses a (D-1)-dimensional pencil decomposition. For 3D:
   * - Real space: distributed in dimensions [0,1], full in dimension [2]
   * - Complex space: full in dimension [0], distributed in dimensions [1,2]
   *
   * Note: We use the FFTW flag for FFTLayoutStruct compatibility since
   * parafaft uses FFTW conventions for r2c padding.
   *
   * Unit test: ctest -R test-parafaftmemorylayout
   **/
  template <size_t NDim> class ParafaftMemoryLayout : public ParafaftPlanner<NDim>
  {
  public:
    ParafaftMemoryLayout() {}

    /**
     * @brief Compute local memory sizes for this MPI rank.
     *
     * Creates a temporary parafaft object to query local sizes, then
     * populates FFTLayoutStruct with the configuration and Fourier space layouts.
     */
    virtual FFTLayoutStruct<NDim> computeLocalSizes(MPICartesianGroup group, device::IdxArray<NDim> nGridPoints,
                                                    bool forbidTransposition = false) override
    {
      // Create FFTLayoutStruct - use FFTW mode for r2c padding compatibility
      // (parafaft uses same padding convention as FFTW)
      FFTLayoutStruct<NDim> result(nGridPoints);

      // Initialize arrays for local layout
      device::IdxArray<NDim> confLocalSizes{};
      device::IdxArray<NDim> confLocalStarts{};
      device::IdxArray<NDim> fourLocalSizes{};
      device::IdxArray<NDim> fourLocalStarts{};
      device::IdxArray<NDim> fourTransposition{};
      device::array<device::IdxArray<2>, NDim> confPadding{};
      std::iota(fourTransposition.begin(), fourTransposition.end(), 0);

      device::Idx parafaftRequiredMemory = 0;

#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
      // Create temporary parafaft object to query sizes
      int globalShape[NDim];
      for (size_t i = 0; i < NDim; ++i)
        globalShape[i] = static_cast<int>(nGridPoints[i]);

      // Use the base communicator - parafaft will create its own Cartesian topology
      parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend> temp(globalShape, group.getBaseComm());

      // Regression guard: if the MPI group was built via FFTMPIDomainSplit, its shape was
      // derived from ParafaftInterface::decomposition — the same probe used here — so these
      // dims must agree. A mismatch means either the user bypassed FFTMPIDomainSplit::makeMPIGroup
      // with a hand-built group (FFTLibrarySelector's verifyDecompositionMatchesBackend should
      // have caught that earlier), or ParaFaFT's probe decision is non-deterministic across
      // planner objects built from the same (baseComm.size(), nGridPoints) — which would be a
      // ParaFaFT-side bug.
      const auto &decomposition = group.getDecomposition();
      int parafaftDecomposition[NDim];
      temp.get_domain_decomposition(parafaftDecomposition);
      for (size_t i = 0; i < NDim; ++i) {
        if (decomposition[i] != parafaftDecomposition[i]) {
          throw ParafaftMemoryLayoutException(
              "ParaFaFT probe disagrees with the MPICartesianGroup shape at dimension ", i, ": probe says ",
              parafaftDecomposition[i], ", group has ", decomposition[i],
              ". Build the group via FFTMPIDomainSplit::makeMPIGroup(baseComm, nGridPoints); "
              "if you did, this indicates ParaFaFT's decomposition heuristic is not deterministic "
              "for these inputs.");
        }
      }

      // Query real (configuration) space layout
      int realShape[NDim], realStart[NDim];
      temp.get_local_real_shape(realShape);
      temp.get_real_global_start(realStart);

      for (size_t i = 0; i < NDim; ++i) {
        confLocalSizes[i] = realShape[i];
        confLocalStarts[i] = realStart[i];
      }

      confPadding[NDim - 1][1] = 2;

      // Query complex (Fourier) space layout
      int complexShape[NDim], complexStart[NDim];
      temp.get_local_complex_shape(complexShape);
      temp.get_complex_global_start(complexStart);

      for (size_t i = 0; i < NDim; ++i) {
        fourLocalSizes[i] = complexShape[i];
        fourLocalStarts[i] = complexStart[i];
      }

      // Memory requirement
      parafaftRequiredMemory = temp.get_required_output_size();
#else
      // Non-parafaft fallback (shouldn't happen)
      for (size_t i = 0; i < NDim; ++i) {
        confLocalSizes[i] = nGridPoints[i];
        fourLocalSizes[i] = nGridPoints[i];
      }
      fourLocalSizes[NDim - 1] = nGridPoints[NDim - 1] / 2 + 1;
      // That's the padding for r2c/cr2, just like in FFTW.
      confPadding[NDim - 1][1] = 2;
#endif
#else
      // Non-MPI fallback (shouldn't happen since parafaft requires MPI)
      for (size_t i = 0; i < NDim; ++i) {
        confLocalSizes[i] = nGridPoints[i];
        fourLocalSizes[i] = nGridPoints[i];
      }
      fourLocalSizes[NDim - 1] = nGridPoints[NDim - 1] / 2 + 1;
      // That's the padding for r2c/cr2, just like in FFTW.
      confPadding[NDim - 1][1] = 2;
#endif

      // Populate result
      result.configurationSpace.setLocalSizes(confLocalSizes);
      result.configurationSpace.setLocalStarts(confLocalStarts);
      result.configurationSpace.setPadding(confPadding);

      result.fourierSpace.setLocalSizes(fourLocalSizes);
      result.fourierSpace.setLocalStarts(fourLocalStarts);
      result.fourierSpace.setTranspositionMap_memoryToGlobalSpace(fourTransposition);

      // Add memory requirement (already in real (double) units)
      result.addExternalMemoryRequest(parafaftRequiredMemory);

      // Set Hermitian partners (same as FFTW)
      result.fourierSpace.setHermitianPartners(
          FFTWHermitianPartners<NDim>::create(result.configurationSpace.getGlobalSizes()));

      return result;
    }
  };
} // namespace TempLat

#endif
