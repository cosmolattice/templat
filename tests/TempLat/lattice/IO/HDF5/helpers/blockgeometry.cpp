
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026

#include <sstream>
#include <string>
#include <vector>

#include "TempLat/lattice/IO/HDF5/helpers/blockgeometry.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct BlockGeometryTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace
  {
    using Idx = device::Idx;

    /** @brief The oracle: the dataset indices the recursive saveDim used to write, reproduced
     *  literally from filesaverhdf5.h before this helper existed.
     *
     *  The non-last dimensions came from the recursion's loop and filter
     *      for (i = 0; i < sizes[dim]; ++i)
     *        if (starts[dim]+i < inicoord || starts[dim]+i >= endcoord ||
     *            (starts[dim]+i - inicoord) % stepcoord != 0) continue;
     *  with the index written being (coord - mdown[dim]) / mstep[dim], i.e. (c-ini)/step.
     *
     *  The last dimension came from the rod branch instead: it is never MPI-split, so the rank
     *  wrote the whole sampled axis, 0 .. (end-ini)/step - 1.
     *
     *  Indices at or beyond the dataset extent are then dropped. Those are the slabs the old code
     *  handed to H5Sselect_hyperslab out of range: the call failed unchecked, the dataspace stayed
     *  at H5S_SELECT_ALL and the H5Dwrite failed too, so nothing of them ever reached the file.
     */
    std::vector<Idx> referenceOffsets(Idx start, Idx size, Idx ini, Idx end, Idx step, bool lastDim, Idx extent)
    {
      std::vector<Idx> raw;
      if (lastDim) {
        const Idx sparseCount = (end - ini) / step;
        for (Idx k = 0; k < sparseCount; ++k)
          raw.push_back(k);
      } else {
        for (Idx i = 0; i < size; ++i) {
          const Idx c = start + i;
          if (c < ini || c >= end || (c - ini) % step != 0) continue;
          raw.push_back((c - ini) / step);
        }
      }

      std::vector<Idx> clamped;
      for (auto k : raw)
        if (k < extent) clamped.push_back(k);
      return clamped;
    }

    std::string describe(const std::string &name, Idx start, Idx size, Idx ini, Idx end, Idx step, Idx extent)
    {
      std::ostringstream oss;
      oss << name << " [start=" << start << " size=" << size << " ini=" << ini << " end=" << end << " step=" << step
          << " extent=" << extent << "]";
      return oss.str();
    }

    /** @brief Check one dimension of the helper against the oracle. */
    void checkDim(TDDAssertion &tdd, const std::string &name, Idx start, Idx size, Idx ini, Idx end, Idx step,
                  Idx extent, bool lastDim = false)
    {
      const auto expected = referenceOffsets(start, size, ini, end, step, lastDim, extent);
      const auto got = computeBlockGeometryDim(start, size, ini, end, step, extent);
      const std::string what = describe(name, start, size, ini, end, step, extent);

      tdd.verify(got.count == (Idx)expected.size(), what + ": count matches the reference sample count");
      if (got.count != (Idx)expected.size()) {
        say << what << ": helper count " << got.count << ", reference count " << expected.size() << "\n";
        return;
      }
      if (expected.empty()) return;

      // The reference list must be the contiguous run [offset, offset+count) the helper claims,
      // which is what makes a single stride-1 hyperslab a faithful replacement for the old
      // one-write-per-rod selection.
      bool sameRun = (got.offset == expected.front());
      for (size_t k = 0; k < expected.size(); ++k)
        sameRun = sameRun && (expected[k] == got.offset + (Idx)k);
      tdd.verify(sameRun, what + ": the reference indices are exactly [offset, offset+count)");

      // And the first coordinate must be the one that lands at that offset.
      tdd.verify(got.first == ini + got.offset * step, what + ": first coordinate matches its dataset offset");
      tdd.verify(got.first >= start && got.first < start + size, what + ": first coordinate is owned by this rank");
    }

    /** @brief Check a whole rank's block, going through the container-level entry point. */
    template <size_t NDim>
    void checkBlock(TDDAssertion &tdd, const std::string &name, const device::IdxArray<NDim> &starts,
                    const device::IdxArray<NDim> &sizes, const device::IdxArray<NDim> &nGrid, bool sparsesave,
                    const std::vector<Idx> &down, const std::vector<Idx> &up, const std::vector<Idx> &step,
                    const std::vector<Idx> &datasetSizes)
    {
      const auto geo = computeSaveBlockGeometry<NDim>(starts, sizes, nGrid, sparsesave, down, up, step, datasetSizes);

      Idx expectedTotal = 1;
      bool expectedEmpty = false;
      for (size_t d = 0; d < NDim; ++d) {
        const Idx ini = sparsesave ? down[d] : 0;
        const Idx end = sparsesave ? up[d] : nGrid[d];
        const Idx stp = sparsesave ? step[d] : 1;
        const bool lastDim = (d == NDim - 1);
        const auto expected = referenceOffsets(starts[d], sizes[d], ini, end, stp, lastDim, datasetSizes[d]);

        std::ostringstream dimName;
        dimName << name << " dim " << d;
        checkDim(tdd, dimName.str(), starts[d], sizes[d], ini, end, stp, datasetSizes[d], lastDim);

        tdd.verify(geo.step[d] == stp, dimName.str() + ": step recorded");
        expectedTotal *= (Idx)expected.size();
        expectedEmpty = expectedEmpty || expected.empty();
      }

      tdd.verify(geo.isEmpty() == expectedEmpty, name + ": emptiness matches the reference");
      tdd.verify(geo.totalCount() == (expectedEmpty ? 0 : expectedTotal), name + ": total element count matches");
    }
  } // namespace

  void BlockGeometryTester::Test(TDDAssertion &tdd)
  {
    const Idx N = 16;

    // ---------------------------------------------------------------------------------------
    // Non-sparse, every MPI split of a 16-point axis. This is the configuration every
    // production checkpoint uses, so it is the one that has to be exactly today's behaviour.
    // ---------------------------------------------------------------------------------------
    for (Idx nRanks : {(Idx)1, (Idx)2, (Idx)4, (Idx)8}) {
      const Idx local = N / nRanks;
      for (Idx r = 0; r < nRanks; ++r) {
        std::ostringstream name;
        name << "non-sparse " << nRanks << "-way rank " << r;
        checkDim(tdd, name.str(), r * local, local, 0, N, 1, N);
      }
    }

    // The last dimension is never MPI-split (ffttopology guarantees it), so it always looks like
    // the one-rank case -- and its offset is 0, which is what the old code hardcoded.
    checkDim(tdd, "non-sparse last dim", 0, N, 0, N, 1, N, /*lastDim=*/true);
    tdd.verify(computeBlockGeometryDim(0, N, 0, N, 1, N).offset == 0, "last dimension starts at file offset 0");

    // ---------------------------------------------------------------------------------------
    // Sparse, step dividing (up-down) exactly: down=2, up=14, step=2 -> extent 6.
    // ---------------------------------------------------------------------------------------
    {
      const Idx down = 2, up = 14, step = 2, extent = (up - down) / step;
      for (Idx nRanks : {(Idx)1, (Idx)2, (Idx)4, (Idx)8}) {
        const Idx local = N / nRanks;
        for (Idx r = 0; r < nRanks; ++r) {
          std::ostringstream name;
          name << "sparse divisible " << nRanks << "-way rank " << r;
          checkDim(tdd, name.str(), r * local, local, down, up, step, extent);
        }
      }
      checkDim(tdd, "sparse divisible last dim", 0, N, down, up, step, extent, /*lastDim=*/true);
    }

    // ---------------------------------------------------------------------------------------
    // Sparse, step NOT dividing (up-down): down=2, up=15, step=2 -> extent 6, but coordinate 14
    // would want index 6. That slab was never written before (out-of-range hyperslab) and is
    // dropped here; the resulting bytes are the same, only the HDF5 error stack disappears.
    // ---------------------------------------------------------------------------------------
    {
      const Idx down = 2, up = 15, step = 2, extent = (up - down) / step;
      tdd.verify(extent == 6, "non-divisible sparse extent floors to 6");
      for (Idx nRanks : {(Idx)1, (Idx)2, (Idx)4, (Idx)8}) {
        const Idx local = N / nRanks;
        for (Idx r = 0; r < nRanks; ++r) {
          std::ostringstream name;
          name << "sparse non-divisible " << nRanks << "-way rank " << r;
          checkDim(tdd, name.str(), r * local, local, down, up, step, extent);
        }
      }
      checkDim(tdd, "sparse non-divisible last dim", 0, N, down, up, step, extent, /*lastDim=*/true);

      // The dropped coordinate is the one the old code could not write: index 6 == extent.
      const auto rank3of4 = computeBlockGeometryDim(12, 4, down, up, step, extent);
      tdd.verify(rank3of4.count == 1, "the rank owning 12..15 keeps only coordinate 12");
      tdd.verify(rank3of4.first == 12 && rank3of4.offset == 5, "coordinate 12 lands at the last valid index");
    }

    // ---------------------------------------------------------------------------------------
    // Ranks that own nothing the dataset wants.
    // ---------------------------------------------------------------------------------------
    {
      // down beyond this rank's range.
      checkDim(tdd, "empty: down > start+size", 0, 4, 8, 12, 1, 4);
      // up below this rank's range.
      checkDim(tdd, "empty: up <= start", 12, 4, 2, 8, 1, 6);
      const auto below = computeBlockGeometryDim(0, 4, 8, 12, 1, 4);
      const auto above = computeBlockGeometryDim(12, 4, 2, 8, 1, 6);
      tdd.verify(below.count == 0 && above.count == 0, "ranks outside the saved window write nothing");
      // A rank straddling the lower limit keeps only the part inside it.
      checkDim(tdd, "partial overlap at the lower limit", 0, 8, 5, 12, 1, 7);
      // A rank straddling the upper limit likewise.
      checkDim(tdd, "partial overlap at the upper limit", 8, 8, 5, 12, 1, 7);
    }

    // ---------------------------------------------------------------------------------------
    // Whole blocks, NDim 1, 2 and 3, through the container entry point.
    // ---------------------------------------------------------------------------------------
    {
      const std::vector<Idx> noDown(10, 0), noUp, noStep(10, 1); // as FileSaverHDF5 default-constructs them
      const std::vector<Idx> full1{N}, full2{N, N}, full3{N, N, N};

      checkBlock<1>(tdd, "NDim 1 non-sparse", {{0}}, {{N}}, {{N}}, false, noDown, noUp, noStep, full1);

      // NDim 2, split 4 ways along dimension 0 only.
      for (Idx r = 0; r < 4; ++r) {
        std::ostringstream name;
        name << "NDim 2 non-sparse rank " << r;
        checkBlock<2>(tdd, name.str(), {{r * 4, 0}}, {{4, N}}, {{N, N}}, false, noDown, noUp, noStep, full2);
      }

      // NDim 3, a pencil decomposition: split along dimensions 0 and 1, never the last.
      for (Idx r0 = 0; r0 < 2; ++r0)
        for (Idx r1 = 0; r1 < 4; ++r1) {
          std::ostringstream name;
          name << "NDim 3 non-sparse rank (" << r0 << "," << r1 << ")";
          checkBlock<3>(tdd, name.str(), {{r0 * 8, r1 * 4, 0}}, {{8, 4, N}}, {{N, N, N}}, false, noDown, noUp, noStep,
                        full3);
        }

      // NDim 3 sparse, the mixed divisible/non-divisible configuration the HDF5 tester saves.
      const std::vector<Idx> down{2, 0, 1}, up{14, N, 15}, step{2, 1, 2};
      std::vector<Idx> sparseSizes(3);
      for (size_t d = 0; d < 3; ++d)
        sparseSizes[d] = (up[d] - down[d]) / step[d];

      for (Idx r0 = 0; r0 < 2; ++r0)
        for (Idx r1 = 0; r1 < 4; ++r1) {
          std::ostringstream name;
          name << "NDim 3 sparse rank (" << r0 << "," << r1 << ")";
          checkBlock<3>(tdd, name.str(), {{r0 * 8, r1 * 4, 0}}, {{8, 4, N}}, {{N, N, N}}, true, down, up, step,
                        sparseSizes);
        }

      // A decomposition where one rank owns nothing of the saved subvolume.
      const std::vector<Idx> tightDown{8, 0, 0}, tightUp{12, N, N}, tightStep{1, 1, 1};
      std::vector<Idx> tightSizes(3);
      for (size_t d = 0; d < 3; ++d)
        tightSizes[d] = (tightUp[d] - tightDown[d]) / tightStep[d];
      checkBlock<3>(tdd, "NDim 3 sparse, rank fully outside", {{0, 0, 0}}, {{4, N, N}}, {{N, N, N}}, true, tightDown,
                    tightUp, tightStep, tightSizes);
      checkBlock<3>(tdd, "NDim 3 sparse, rank fully inside", {{8, 0, 0}}, {{4, N, N}}, {{N, N, N}}, true, tightDown,
                    tightUp, tightStep, tightSizes);
    }

    // ---------------------------------------------------------------------------------------
    // The loader's configuration: the dense one. FileLoaderHDF5 has no sparse support, so its
    // dataset covers the lattice one to one and the arithmetic must collapse to
    // offset = first = starts[d], count = sizes[d] -- the raw global coordinates the recursive
    // loadDim used as offsets.
    // ---------------------------------------------------------------------------------------
    {
      const device::IdxArray<3> nGrid{{16, 16, 16}};
      for (Idx r0 = 0; r0 < 2; ++r0)
        for (Idx r1 = 0; r1 < 4; ++r1) {
          const device::IdxArray<3> starts{{r0 * 8, r1 * 4, 0}}, sizes{{8, 4, 16}};
          const auto geo = computeDenseBlockGeometry<3>(starts, sizes, nGrid);
          bool ok = true;
          for (size_t d = 0; d < 3; ++d)
            ok = ok && geo.count[d] == sizes[d] && geo.first[d] == starts[d] && geo.offset[d] == starts[d] &&
                 geo.step[d] == 1;
          std::ostringstream name;
          name << "dense (loader) rank (" << r0 << "," << r1 << ")";
          tdd.verify(ok, name.str() + ": covers exactly this rank's subdomain at its global offset");
          tdd.verify(geo.totalCount() == 8 * 4 * 16, name.str() + ": stages the whole local block");
        }

      // Non-cubic, to catch a dimension swap in the loader too.
      const device::IdxArray<3> odd{{8, 16, 32}};
      const auto geo = computeDenseBlockGeometry<3>(device::IdxArray<3>{{4, 8, 0}}, device::IdxArray<3>{{4, 8, 32}},
                                                    odd);
      tdd.verify(geo.count[0] == 4 && geo.count[1] == 8 && geo.count[2] == 32,
                 "dense (loader) non-cubic: extents follow the lattice dimension order");
      tdd.verify(geo.offset[0] == 4 && geo.offset[1] == 8 && geo.offset[2] == 0,
                 "dense (loader) non-cubic: offsets are the global starts, last dimension at 0");
    }

    // ---------------------------------------------------------------------------------------
    // Non-cubic lattice: a dimension swap or a shared stride would show up here and nowhere else.
    // ---------------------------------------------------------------------------------------
    {
      const std::vector<Idx> noDown(10, 0), noUp, noStep(10, 1);
      const std::vector<Idx> sizes{8, 16, 32};
      for (Idx r = 0; r < 4; ++r) {
        std::ostringstream name;
        name << "non-cubic 8x16x32 rank " << r;
        checkBlock<3>(tdd, name.str(), {{r * 2, 0, 0}}, {{2, 16, 32}}, {{8, 16, 32}}, false, noDown, noUp, noStep,
                      sizes);
      }
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::BlockGeometryTester> test;
}
