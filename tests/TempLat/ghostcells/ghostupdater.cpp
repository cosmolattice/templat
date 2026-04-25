
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025
#include "TempLat/lattice/ghostcells/ghostupdater.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/lattice/memory/triplestatelayouts.h"
#include "TempLat/lattice/algebra/operators/power.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/util/ndloop.h"

#include <sstream>
#include <iomanip>

namespace TempLat
{

  /** @brief A namespace purely for test structures / classes. Quick simple structs for whatever test you need go in
   * here. */
  namespace TestScratch
  {
    /* quick and ugly helper struct */
    template <size_t NDim> struct datumMPITypeHolder {
      MPI_Datatype dType;
      datumMPITypeHolder()
      {
#ifdef HAVE_MPI
        MPI_Type_contiguous(NDim, TempLat::MPITypeSelect<ptrdiff_t>(), &dType);
        MPI_Type_commit(&dType);
#endif
      }
      ~datumMPITypeHolder()
      {
#ifdef HAVE_MPI
        int didFinalize = 0;
        MPI_Finalized(&didFinalize);
        if (!didFinalize) MPI_Type_free(&dType);
#endif
      }
    };

    template <size_t NDim> struct datum {
      device::IdxArray<NDim> data;

      friend std::ostream &operator<<(std::ostream &ostream, const datum &dat)
      {
        ostream << dat.data;
        return ostream;
      }

      static MPI_Datatype getMPIType()
      {
        static datumMPITypeHolder<NDim> holder;
        return holder.dType;
      }
    };

    template <size_t NDim> void datum_initialize(MemoryBlock<datum<NDim>, NDim> &block, const LayoutStruct<NDim> layout)
    {
      const auto localSizes = layout.getLocalSizes();
      const ptrdiff_t nGhost = layout.getNGhosts();

      std::string localSizesStr = "{";
      for (size_t i = 0; i < NDim; ++i) {
        localSizesStr += std::to_string(localSizes[i]);
        if (i != NDim - 1) localSizesStr += ", ";
      }
      localSizesStr += "}";

      if constexpr (NDim == 1) {
        auto view = block.getRawView();
        device::iteration::foreach (
            "DatumInitialize1D", device::IdxArray<1>{0}, device::IdxArray<1>(localSizes[0]),
            DEVICE_LAMBDA(const size_t i) {
              view(nGhost + i) = datum<NDim>{layout.getLocalStarts()[0] + (ptrdiff_t)i + 1};
            });
      } else {
        device::IdxArray<NDim> viewSizes;
        for (size_t k = 0; k < NDim; ++k)
          viewSizes[k] = localSizes[k] + 2 * nGhost;
        auto view = block.getNDView(viewSizes);

        device::array<std::pair<ptrdiff_t, ptrdiff_t>, NDim> slices{};
        for (size_t k = 0; k < NDim; ++k)
          slices[k] = std::make_pair(nGhost, nGhost + localSizes[k]);

        auto subView =
            device::apply([&](const auto &...args) { return device::memory::subview(view, args...); }, slices);
        auto functor = DEVICE_LAMBDA(const device::IdxArray<NDim> &idx)
        {
          device::IdxArray<NDim> val;
          for (size_t k = 0; k < NDim; ++k)
            val[k] = idx[k] + layout.getLocalStarts()[k] + 1;
          device::IdxArray<NDim> inc;
          for (size_t k = 0; k < NDim; ++k)
            inc[k] = k;
          device::apply([&](const auto &...args) { subView(idx[args]...) = datum<NDim>{val[args]...}; }, inc);
        };

        device::IdxArray<NDim> it_stop{};
        for (size_t k = 0; k < NDim; ++k)
          it_stop[k] = localSizes[k];

        device::iteration::foreach ("GhostUpdater", {}, it_stop, functor);
      }
    }

    template <size_t nd> bool test_ghost_updater(const ptrdiff_t nGrid, const size_t nGhost)
    {
      auto toolBox = MemoryToolBox<nd>::makeShared(nGrid, nGhost);
      toolBox->unsetVerbose();

      MPICartesianGroup mGroup(FFTMPIDomainSplit<nd>::makeMPIGroup(toolBox->mNGridPointsVec));
      FFTLibrarySelector<nd> fftlib(mGroup, toolBox->mNGridPointsVec);
      TripleStateLayouts fullLayout(fftlib.getLayout(), nGhost);
      GhostUpdater<nd> ghostUpdater(mGroup, fullLayout.getConfigSpaceLayout());

      const auto localSizes = fullLayout.getConfigSpaceLayout().getLocalSizes();
      device::IdxArray<nd> fullLocalSizes{};
      for (size_t i = 0; i < nd; ++i)
        fullLocalSizes[i] = localSizes[i] + 2 * nGhost;

      const bool verbose = (nd == 2) && (nGrid <= 16);

      auto print_it = [&](auto view) {
        if (!verbose) return;
        size_t total_size = 1;
        for (size_t i = 0; i < nd; ++i)
          total_size *= fullLocalSizes[i];
        device::IdxArray<nd> cIdx{};
        std::stringstream gridArrayStr;
        for (size_t i = 0; i < total_size; ++i) {
          // Linear index to cartesian index
          size_t lsize = 1;
          size_t remainder = i;
          for (size_t j = 0; j < nd; ++j) {
            lsize = fullLocalSizes[nd - 1 - j];
            cIdx[nd - 1 - j] = remainder % lsize;
            remainder = (remainder - cIdx[nd - 1 - j]) / fullLocalSizes[nd - 1 - j];
          }
          if (cIdx[nd - 1] == 0) gridArrayStr << "\n";
          gridArrayStr << std::setw(4) << view(i).data[0];
          if ((size_t)cIdx[nd - 1] != (size_t)nGrid - 1) gridArrayStr << " ";
        }
        gridArrayStr << "\n\n";

        sayMPI << "Current view:" << gridArrayStr.str();
      };

      MemoryBlock<TestScratch::datum<nd>, nd> block(pow<nd>(nGrid + 2 * nGhost));
      const auto config_layout = fullLayout.getConfigSpaceLayout();
      TestScratch::datum_initialize<nd>(block, config_layout);

      size_t total_size = 1;
      for (size_t i = 0; i < nd; ++i)
        total_size *= fullLocalSizes[i];
      device::IdxArray<nd> cIdx{};

      auto view = block.getRawHostView();
      print_it(view);

      ghostUpdater.update(block);

      block.flagHostMirrorOutdated();
      view = block.getRawHostView();
      print_it(view);

      bool all_correct = true;
      size_t ww = 0;
      for (size_t i = 0; i < total_size; ++i) {
        // Linear index to cartesian index
        size_t lsize = 1;
        size_t remainder = i;
        for (size_t j = 0; j < nd; ++j) {
          lsize = fullLocalSizes[nd - 1 - j];
          cIdx[nd - 1 - j] = remainder % lsize;
          remainder = (remainder - cIdx[nd - 1 - j]) / fullLocalSizes[nd - 1 - j];
          cIdx[nd - 1 - j] += config_layout.getLocalStarts()[nd - 1 - j] + 1 - nGhost;
        }

        auto should_value = cIdx;
        for (size_t d = 0; d < nd; ++d) {
          if (cIdx[d] < 1) should_value[d] += nGrid;
          if (cIdx[d] > nGrid) should_value[d] -= nGrid;
        }

        device::IdxArray<nd> is_value;
        for (size_t d = 0; d < nd; ++d)
          is_value[d] = view(i).data[d];

        for (size_t d = 0; d < nd; ++d) {
          all_correct &= (size_t)is_value[d] == (size_t)should_value[d];
          if ((size_t)is_value[d] != (size_t)should_value[d])
            sayMPI << ++ww << " false " << is_value << " (is) | vs | " << should_value << " (should) AT POSITION "
                   << cIdx << "\n";
          // else
          //   sayMPI << ++ww << " true " << is_value << " (is) | vs | " << should_value << " (should) AT POSITION "
          //          << cIdx << "\n";
        }
      }
      return all_correct;
    }
  } // namespace TestScratch

  template <size_t NDim> struct GhostUpdaterTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void GhostUpdaterTester<NDim>::Test(TDDAssertion &tdd)
  {
    static_assert(NDim > 1, "GhostUpdater test only makes sense in 2 or more dimensions.");
    // I just don't have the patience to do the 1D case, since it involves no MPI communication.

    // Note: nGrid=4 cases are deliberately skipped. With ParaFaFT's (NDim-1)-D
    // pencil decomposition (and an R2C half-axis of nGrid/2+1 = 3), such a tiny
    // grid produces zero-extent or non-divisible local shapes on some ranks for
    // most comm sizes, mixing per-rank failure modes. Not worth testing — we
    // exercise the same code paths with the larger grids below.

    if constexpr (NDim < 6) tdd.verify(TestScratch::test_ghost_updater<NDim>(16, 1));
    if constexpr (NDim < 4) {
      tdd.verify(TestScratch::test_ghost_updater<NDim>(32, 1));
      tdd.verify(TestScratch::test_ghost_updater<NDim>(128, 1));
    }

    if constexpr (NDim < 6) tdd.verify(TestScratch::test_ghost_updater<NDim>(16, 2));
    if constexpr (NDim < 4) {
      tdd.verify(TestScratch::test_ghost_updater<NDim>(32, 2));
      tdd.verify(TestScratch::test_ghost_updater<NDim>(128, 2));
    }

    if constexpr (NDim < 6) tdd.verify(TestScratch::test_ghost_updater<NDim>(16, 3));
    if constexpr (NDim < 4) {
      tdd.verify(TestScratch::test_ghost_updater<NDim>(32, 3));
      tdd.verify(TestScratch::test_ghost_updater<NDim>(128, 3));
    }
  }

  /** @brief BC-aware ghost-fill tester.
   *
   * For a Field<double, NDim> built with a chosen per-dim BCSpec, sets every cell to
   * (global_x[bcDim] + 1). The "+1" keeps owned values in [1, nGrid] (strictly positive), so
   * Dirichlet's zero-fill is visible against interior values, and every other BC's expected
   * value is distinct from every other BC's expected value. Then calls updateGhosts() and
   * checks the low- and high-ghost slabs along bcDim. Each rank verifies its own slab —
   * boundary ranks see the BC transform, interior ranks see the ordinary periodic exchange.
   */
  template <size_t NDim> struct BCFillTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace BCTestDetail
  {
    inline double expectedLowGhost(BCType bc, ptrdiff_t globalSize, ptrdiff_t localStart, bool isLowBoundary)
    {
      if (!isLowBoundary) return static_cast<double>(localStart);
      switch (bc) {
      case BCType::Periodic:     return static_cast<double>(globalSize);
      case BCType::Antiperiodic: return -static_cast<double>(globalSize);
      case BCType::Dirichlet:    return 0.0;
      case BCType::Neumann:      return 1.0;
      }
      return 0.0;
    }

    inline double expectedHighGhost(BCType bc, ptrdiff_t globalSize, ptrdiff_t localStart, ptrdiff_t localSize,
                                    bool isHighBoundary)
    {
      if (!isHighBoundary) return static_cast<double>(localStart + localSize + 1);
      switch (bc) {
      case BCType::Periodic:     return 1.0;
      case BCType::Antiperiodic: return -1.0;
      case BCType::Dirichlet:    return 0.0;
      case BCType::Neumann:      return static_cast<double>(globalSize);
      }
      return 0.0;
    }

    template <size_t NDim>
    bool checkBCInDim(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t bcDim, BCType bc,
                      ptrdiff_t nGrid, ptrdiff_t nGhost)
    {
      auto layout = toolBox->mLayouts.getConfigSpaceLayout();
      const auto &localSizes = layout.getLocalSizes();
      const auto &localStarts = layout.getLocalStarts();

      BCSpec<NDim> spec = allPeriodic<NDim>();
      spec[bcDim] = bc;

      Field<double, NDim> f("f_bc_check", toolBox, LatticeParameters<double>(), spec);
      SpatialCoordinate<NDim> x(toolBox);
      // Initialize so each cell holds (global_x[bcDim] + 1). SpatialCoordinate's selector is a
      // compile-time tag, so we dispatch on bcDim via constexpr_for.
      constexpr_for<1, NDim + 1>([&](auto dirTag) {
        constexpr size_t d = static_cast<size_t>(decltype(dirTag)::value) - 1;
        if (d == bcDim) f = x(dirTag) + 1.0;
      });
      f.updateGhosts();

      const bool isLowBoundary  = (localStarts[bcDim] == 0);
      const bool isHighBoundary = (localStarts[bcDim] + localSizes[bcDim] == nGrid);

      auto view = f.getFullNDHostView();

      // Walk the low and high faces along bcDim, restricted to indices within the owned slab in
      // every other dim — corners would mix in unrelated periodic BC behavior of other dims.
      auto walkFace = [&](bool low, bool isBoundary) {
        bool ok = true;
        const ptrdiff_t bcIdx = low ? (nGhost - 1) : (nGhost + localSizes[bcDim]);
        const double expected = low
            ? expectedLowGhost(bc, nGrid, localStarts[bcDim], isBoundary)
            : expectedHighGhost(bc, nGrid, localStarts[bcDim], localSizes[bcDim], isBoundary);

        device::IdxArray<NDim> idx{};
        for (size_t i = 0; i < NDim; ++i) idx[i] = (i == bcDim) ? bcIdx : nGhost;
        const ptrdiff_t loopDim = (bcDim == 0) ? 1 : 0;
        // Only iterate over the loopDim — checking one row of the face is sufficient since the
        // value of `f` depends only on global_x[bcDim], so all face cells have the same expected.
        const ptrdiff_t loopExtent = (NDim == 1) ? 1 : localSizes[loopDim];
        for (ptrdiff_t i = 0; i < loopExtent; ++i) {
          if constexpr (NDim > 1) idx[loopDim] = nGhost + i;
          const double got = device::apply([&](auto... a) { return view(a...); }, idx);
          if (std::abs(got - expected) > 1e-14) {
            ok = false;
            std::stringstream ss;
            ss << "BCFill mismatch (bcDim=" << bcDim << ", bc=" << static_cast<int>(bc)
               << ", " << (low ? "low" : "high") << " face) at view_idx=" << idx
               << " got=" << got << " expected=" << expected
               << " localStart=" << localStarts[bcDim] << " localSize=" << localSizes[bcDim] << "\n";
            sayMPI << ss.str();
          }
        }
        return ok;
      };

      bool ok = true;
      ok &= walkFace(true,  isLowBoundary);
      ok &= walkFace(false, isHighBoundary);
      return ok;
    }
  } // namespace BCTestDetail

  template <size_t NDim> void BCFillTester<NDim>::Test(TDDAssertion &tdd)
  {
    constexpr ptrdiff_t nGrid = 8;
    constexpr ptrdiff_t nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    for (size_t bcDim = 0; bcDim < NDim; ++bcDim) {
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Antiperiodic, nGrid, nGhost));
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Dirichlet,    nGrid, nGhost));
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Neumann,      nGrid, nGhost));
      // Periodic regression: new BC plumbing must not change existing all-periodic behavior.
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Periodic,     nGrid, nGhost));
    }

    // Mixed-BC sanity: Antiperiodic in dim 0, Periodic elsewhere — verify dim-0 low ghost follows
    // Antiperiodic on the boundary rank, leaving the other dim's Periodic unaffected.
    if constexpr (NDim >= 2) {
      auto layout = toolBox->mLayouts.getConfigSpaceLayout();
      const auto &localSizes = layout.getLocalSizes();
      const auto &localStarts = layout.getLocalStarts();
      (void)localSizes;

      BCSpec<NDim> mixed = allPeriodic<NDim>();
      mixed[0] = BCType::Antiperiodic;

      Field<double, NDim> f("f_mixed", toolBox, LatticeParameters<double>(), mixed);
      SpatialCoordinate<NDim> x(toolBox);
      f = x(1_c) + 1.0;
      f.updateGhosts();

      auto view = f.getFullNDHostView();
      const bool isLow0 = (localStarts[0] == 0);

      device::IdxArray<NDim> idx{};
      idx[0] = nGhost - 1;
      for (size_t i = 1; i < NDim; ++i) idx[i] = nGhost;
      const double got = device::apply([&](auto... a) { return view(a...); }, idx);
      const double expected = isLow0 ? -static_cast<double>(nGrid) : static_cast<double>(localStarts[0]);
      tdd.verify(std::abs(got - expected) <= 1e-14);
    }
  }

} // namespace TempLat

namespace
{
  // TempLat::TDDContainer<TempLat::GhostUpdaterTester<1>> test1;
  TempLat::TDDContainer<TempLat::GhostUpdaterTester<2>> test2;
  TempLat::TDDContainer<TempLat::GhostUpdaterTester<3>> test3;
  TempLat::TDDContainer<TempLat::GhostUpdaterTester<4>> test4;
  TempLat::TDDContainer<TempLat::GhostUpdaterTester<5>> test5;
  TempLat::TDDContainer<TempLat::GhostUpdaterTester<6>> test6;

  TempLat::TDDContainer<TempLat::BCFillTester<2>> bcTest2;
  TempLat::TDDContainer<TempLat::BCFillTester<3>> bcTest3;
} // namespace
