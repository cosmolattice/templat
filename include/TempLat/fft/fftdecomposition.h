#ifndef TEMPLAT_FFT_FFTDECOMPOSITION_H
#define TEMPLAT_FFT_FFTDECOMPOSITION_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

#include <array>
#include <cstddef>

namespace TempLat
{
  /** @brief The MPI Cartesian decomposition an FFT backend will use for a given runtime setup.
   *
   *  Backends advertise this via a static `decomposition(baseComm, nGridPoints)` method so the
   *  MPI Cartesian group can be built to match the backend's requirements before any plan is
   *  constructed.
   *
   *  Fields:
   *  - nDimsToSplit: number of leading dimensions MPI should distribute over. A value of 0 means
   *    the backend does not use an MPI decomposition (e.g. KokkosFFT, single-rank FFTW).
   *  - dims: optional explicit grid shape (rank counts per lattice dimension). An entry of 0 means
   *    "let `MPI_Dims_create` choose". If all nDimsToSplit leading entries are nonzero, the shape
   *    is pinned exactly and the product must equal the rank count.
   */
  template <size_t NDim> struct FFTDecomposition {
    ptrdiff_t nDimsToSplit = 0;
    std::array<int, NDim> dims = {};
  };
} // namespace TempLat

#endif
