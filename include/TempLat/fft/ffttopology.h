#ifndef TEMPLAT_FFT_FFTTOPOLOGY_H
#define TEMPLAT_FFT_FFTTOPOLOGY_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"
#include "TempLat/util/exception.h"

#include <array>
#include <cstddef>
#include <vector>

namespace TempLat
{
  MakeException(FFTTopologyException);

  /** @brief The complete MPI topology an FFT backend will use: not just the grid shape, but the
   *  communicator whose rank order defines it and this rank's position within it.
   *
   *  This exists because shape agreement is not enough. TempLat and the FFT backend used to build
   *  Cartesian communicators independently, both with `reorder = 1`, and reconcile them by
   *  comparing grid *shape* only. But the local starts installed into the layout are derived from
   *  the backend's coordinates, while ghost neighbours come from `MPI_Cart_shift` on TempLat's
   *  communicator. Nothing forced the two rank->coordinate mappings to agree; it worked only
   *  because mainstream MPI implementations ignore `reorder`. When they disagree, subdomain
   *  boundaries silently exchange the wrong data.
   *
   *  So the backend is the single source of truth. It reports the communicator it will actually
   *  plan on, and `MPICartesianGroup` builds its grid over that communicator with
   *  `reorder = false`, then verifies the coordinates it computes match `coords` below.
   *
   *  Fields:
   *  - comm: the communicator the backend plans on. Its rank order is authoritative. For ParaFaFT
   *    this is a duplicate of its internal Cartesian communicator; for FFTW and KokkosFFT, which
   *    build no Cartesian topology of their own, it is the base communicator.
   *  - dims: grid shape, one entry per lattice dimension, trailing 1s for undistributed axes.
   *  - coords: this rank's coordinates in that grid. Undistributed axes are 0.
   *  - nDimsToSplit: number of dimensions actually distributed (entries of `dims` greater than 1).
   */
  template <size_t NDim> struct FFTTopology {
    MPICommReference comm;
    std::array<int, NDim> dims = {};
    std::array<int, NDim> coords = {};
    int nDimsToSplit = 0;

    /** @brief Shape-only view, for the call sites that genuinely only need the grid. */
    FFTDecomposition<NDim> decomposition() const { return FFTDecomposition<NDim>{nDimsToSplit, dims}; }

    std::vector<int> dimsVector() const { return std::vector<int>(dims.begin(), dims.end()); }

    /** @brief Reject a topology that would violate the invariants the rest of TempLat relies on.
     *
     *  Call this at every point a topology is produced, so a backend change fails loudly at
     *  startup rather than corrupting data or output later.
     *
     *  Not named `verify`: that is a function-like macro in TempLat/util/tdd/tddmacros.h.
     */
    void checkInvariants(int commSize) const
    {
      long product = 1;
      for (size_t i = 0; i < NDim; ++i) {
        if (dims[i] < 1)
          throw FFTTopologyException("FFT backend reported a non-positive grid extent ", dims[i], " at dimension ", i,
                                     ".");
        if (coords[i] < 0 || coords[i] >= dims[i])
          throw FFTTopologyException("FFT backend reported coordinate ", coords[i], " at dimension ", i,
                                     ", outside the grid extent ", dims[i], ".");
        product *= dims[i];
      }

      if (product != commSize)
        throw FFTTopologyException("FFT backend grid has ", product, " cells but the communicator has ", commSize,
                                   " ranks. The grid must tile the communicator exactly.");

      // The last axis must stay local. Three places depend on this and none of them can detect a
      // violation on their own:
      //   - ParafaftMemoryLayout applies the r2c padding to the last axis only;
      //   - ParaFaFT pins real_global_start_[D-1] = 0;
      //   - FileSaverHDF5 writes the last axis as one contiguous rod with offset 0, so a split
      //     there corrupts output files with no error at all.
      if (dims[NDim - 1] != 1)
        throw FFTTopologyException(
            "FFT backend wants to distribute the last lattice dimension (grid extent ", dims[NDim - 1],
            "), which TempLat does not support. The r2c padding in ParafaftMemoryLayout, ParaFaFT's "
            "real_global_start_[D-1] = 0, and the contiguous last-axis write in FileSaverHDF5 all assume it "
            "stays local. Fix all three before allowing this.");
    }
  };
} // namespace TempLat

#endif
