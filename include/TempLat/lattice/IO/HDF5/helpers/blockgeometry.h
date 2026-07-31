#ifndef TEMPLAT_LATTICE_IO_HDF5_HELPERS_BLOCKGEOMETRY_H
#define TEMPLAT_LATTICE_IO_HDF5_HELPERS_BLOCKGEOMETRY_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026

#include <cstddef>

#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief What one rank contributes to one dataset dimension.
   *
   *  The dataset stores the sampled coordinates ini, ini+step, ini+2*step, ... at consecutive
   *  indices 0, 1, 2, ...; this rank owns the contiguous index run [offset, offset+count), whose
   *  first element is the global lattice coordinate `first`.
   */
  struct BlockGeometryDim {
    device::Idx count = 0;  ///< number of samples this rank writes along this dimension
    device::Idx first = 0;  ///< first global lattice coordinate this rank writes
    device::Idx offset = 0; ///< index of `first` inside the dataset
  };

  /** @brief The per-dimension hyperslab arithmetic behind FileSaverHDF5::saveBlock.
   *
   *  @param start          global coordinate of this rank's first owned point (getConfigSpaceStarts)
   *  @param size           number of points this rank owns (getConfigSpaceSizes)
   *  @param ini            first coordinate the dataset samples (0, or `down` when saving sparsely)
   *  @param end            one past the last coordinate the dataset samples (nGrid, or `up`)
   *  @param step           sampling stride (1, or `step`); values below 1 are treated as 1
   *  @param datasetExtent  the dataset's extent in this dimension, i.e. (up-down)/step (or nGrid)
   *
   *  The saved coordinates are ini + k*step for k = 0, 1, ...; their file index is exactly k
   *  (that is (c-down)/step, the index the recursive saveDim used to compute point by point), so
   *  the coordinates one rank owns map to *consecutive* file indices and its region is a single
   *  stride-1, block-1 hyperslab.
   *
   *  Clamping to `datasetExtent` matters when step does not divide (end-ini): the last surviving
   *  coordinate would land one past the extent. The old recursion let that slab fail inside HDF5
   *  (an unchecked H5Sselect_hyperslab, hence an unchecked H5Dwrite), so it was never written;
   *  dropping it here silently produces the same bytes without the error stack.
   *
   *  Clipping `hi` at start+size matters when up > nGrid: the old code walked off the end of the
   *  rank's allocation instead. FileSaverHDF5::save rejects such limits outright, so the clip is
   *  belt and braces.
   */
  inline BlockGeometryDim computeBlockGeometryDim(device::Idx start, device::Idx size, device::Idx ini,
                                                  device::Idx end, device::Idx step, device::Idx datasetExtent)
  {
    const device::Idx stp = step < 1 ? 1 : step;

    const device::Idx lo = start > ini ? start : ini;                          // max(start, ini)
    const device::Idx hi = (start + size) < end ? (start + size) : end;        // min(start+size, end)

    // Ceiling division; both numerators are non-negative where they are used.
    const device::Idx kFirst = (lo <= ini) ? 0 : (lo - ini + stp - 1) / stp;
    const device::Idx kEnd = (hi <= ini) ? 0 : (hi - ini + stp - 1) / stp;
    const device::Idx kStop = kEnd < datasetExtent ? kEnd : datasetExtent;

    BlockGeometryDim result;
    result.count = kStop > kFirst ? kStop - kFirst : 0;
    result.first = ini + kFirst * stp;
    result.offset = kFirst;
    return result;
  }

  /** @brief The block of a lattice dataset that the calling rank writes: one hyperslab, whatever
   *  the MPI decomposition and whatever the sparse limits.
   *
   *  Dimensions are in lattice order, which is also file order and memory order: the configuration
   *  space layout is never transposed (only the Fourier layout sets a transposition map), so the
   *  d-th dataset dimension is the d-th memory dimension.
   */
  template <size_t NDim> struct SaveBlockGeometry {
    device::IdxArray<NDim> count{};  ///< extent of this rank's hyperslab, per dimension
    device::IdxArray<NDim> first{};  ///< first global lattice coordinate written, per dimension
    device::IdxArray<NDim> offset{}; ///< offset of the hyperslab in the dataset, per dimension
    device::IdxArray<NDim> step{};   ///< coordinate stride between consecutive samples

    /** @brief True when this rank writes nothing at all, which is legal: the transfer mode is
     *  H5FD_MPIO_INDEPENDENT, so a rank may simply skip its H5Dwrite. */
    bool isEmpty() const
    {
      for (size_t d = 0; d < NDim; ++d)
        if (count[d] <= 0) return true;
      return false;
    }

    device::Idx totalCount() const
    {
      if (isEmpty()) return 0;
      device::Idx total = 1;
      for (size_t d = 0; d < NDim; ++d)
        total *= count[d];
      return total;
    }
  };

  /** @brief Assemble the per-dimension results into one block.
   *
   *  @param starts,sizes  this rank's local decomposition, from TripleStateLayouts
   *  @param nGrid         global lattice size, used as `end` when not saving sparsely
   *  @param sparsesave    whether setLimits was called
   *  @param down,up,step  the sparse limits; `up` is only read when sparsesave is true, as
   *                       FileSaverHDF5 leaves it empty otherwise
   *  @param datasetSizes  the extents the dataset was created with
   */
  template <size_t NDim, typename CStarts, typename CSizes, typename CGrid, typename CLimits>
  SaveBlockGeometry<NDim> computeSaveBlockGeometry(const CStarts &starts, const CSizes &sizes, const CGrid &nGrid,
                                                   bool sparsesave, const CLimits &down, const CLimits &up,
                                                   const CLimits &step, const CLimits &datasetSizes)
  {
    SaveBlockGeometry<NDim> geo{};
    for (size_t d = 0; d < NDim; ++d) {
      const device::Idx ini = sparsesave ? (device::Idx)down[d] : 0;
      const device::Idx end = sparsesave ? (device::Idx)up[d] : (device::Idx)nGrid[d];
      const device::Idx stp = sparsesave ? (device::Idx)step[d] : 1;

      const auto dim = computeBlockGeometryDim((device::Idx)starts[d], (device::Idx)sizes[d], ini, end, stp,
                                               (device::Idx)datasetSizes[d]);
      geo.count[d] = dim.count;
      geo.first[d] = dim.first;
      geo.offset[d] = dim.offset;
      geo.step[d] = stp;
    }
    return geo;
  }

  /** @brief The dense configuration of computeSaveBlockGeometry: the whole lattice, every point.
   *
   *  This is what FileLoaderHDF5 needs -- a dataset that covers the lattice one-to-one, which is
   *  the only thing it can read, since it has no sparse support. With ini = 0, end = nGrid and
   *  step = 1 the arithmetic collapses to offset = first = starts[d] and count = sizes[d].
   */
  template <size_t NDim, typename CStarts, typename CSizes, typename CGrid>
  SaveBlockGeometry<NDim> computeDenseBlockGeometry(const CStarts &starts, const CSizes &sizes, const CGrid &nGrid)
  {
    // The limit containers are never read when sparsesave is false; nGrid stands in so that all
    // four share a type.
    return computeSaveBlockGeometry<NDim>(starts, sizes, nGrid, false, nGrid, nGrid, nGrid, nGrid);
  }

} // namespace TempLat

#endif
