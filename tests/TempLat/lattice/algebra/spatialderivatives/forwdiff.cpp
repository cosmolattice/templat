
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/spatialderivatives/forwdiff.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/helpers/getvectorcomponent.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/util/ndloop.h"
#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{

  template <size_t NDim> struct ForwDiffTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> inline void ForwDiffTester<NDim>::Test(TDDAssertion &tdd)
  {
    const device::Idx nGrid = 8, nGhost = 1;

    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    SpatialCoordinate<NDim> x(toolBox);
    toolBox->setVerbose();

    // Get layout for computing global coordinates (for error reporting)
    auto layout = toolBox->mLayouts.getConfigSpaceLayout();

    // Test forward difference in each dimension
    // ForwDiff<dir> computes: (f[idx + e_dir] - f[idx]) / dx
    constexpr_for<1, NDim + 1>([&](auto dirTag) {
      constexpr int dir = decltype(dirTag)::value;
      constexpr size_t d = static_cast<size_t>(dir) - 1;

      Field<double, NDim> sc("SC_" + std::to_string(d), toolBox);
      sc = x(dirTag);
      sc.updateGhosts();

      Field<double, NDim> fdsc("fdSC_" + std::to_string(d), toolBox);
      fdsc = forwDiff(sc, Tag<dir>{});

      bool OK = true;
      auto sc_view = sc.getFullNDHostView();
      auto fdsc_view = fdsc.getLocalNDHostView();

      // Use NDLoop to iterate
      // For forward diff in dir: fdsc[idx] = (sc[idx+nGhost+e_dir] - sc[idx+nGhost]) / dx
      NDLoop<NDim>(fdsc_view, [&](const auto &...indices) {
        // Build index arrays for accessing full view (with ghosts)
        device::IdxArray<NDim> idx_base = {(indices + nGhost)...};
        device::IdxArray<NDim> idx_next = idx_base;
        idx_next[d] += 1;

        // Get values
        const double val_current = device::apply([&](auto... i) { return sc_view(i...); }, idx_base);
        const double val_next = device::apply([&](auto... i) { return sc_view(i...); }, idx_next);
        const double val_fdsc = fdsc_view(indices...);

        const double dx = 1.0; // dx = 1 by default in these tests
        const double expect = (val_next - val_current) / dx;

        if (std::abs(expect - val_fdsc) > 1e-14) {
          OK = false;

          // Compute global coordinates for error message
          device::IdxArray<NDim> global_idx;
          layout.putSpatialLocationFromMemoryIndexInto(global_idx, indices...);

          std::stringstream ss;
          ss << "ForwDiff mismatch in dir " << d << " at global (";
          for (size_t i = 0; i < NDim; ++i) {
            ss << global_idx[i];
            if (i < NDim - 1) ss << ", ";
          }
          ss << "): expect = " << expect << ", fdSC = " << val_fdsc << "\n";
          sayMPI << ss.str();
        }
      });
      tdd.verify(OK);
    });

    // Test 2: ForwDiff of quadratic function sc^2
    // For f(x) = x^2, forward difference is ((x+1)^2 - x^2) = 2x + 1
    {
      Field<double, NDim> sc("SC_sq", toolBox);
      sc = x(1_c);
      sc.updateGhosts();

      Field<double, NDim> sc_sq("SC_sq_field", toolBox);
      sc_sq = sc * sc;
      sc_sq.updateGhosts();

      Field<double, NDim> fdsc_sq("fdSC_sq", toolBox);
      fdsc_sq = forwDiff(sc * sc, Tag<1>{});

      bool OK = true;
      auto sc_sq_view = sc_sq.getFullNDHostView();
      auto fdsc_sq_view = fdsc_sq.getLocalNDHostView();

      constexpr size_t d = 0;

      NDLoop<NDim>(fdsc_sq_view, [&](const auto &...indices) {
        device::IdxArray<NDim> idx_base = {(indices + nGhost)...};
        device::IdxArray<NDim> idx_next = idx_base;
        idx_next[d] += 1;

        const double val_current = device::apply([&](auto... i) { return sc_sq_view(i...); }, idx_base);
        const double val_next = device::apply([&](auto... i) { return sc_sq_view(i...); }, idx_next);
        const double val_fdsc_sq = fdsc_sq_view(indices...);

        const double dx = 1.0;
        const double expect = (val_next - val_current) / dx;

        if (std::abs(expect - val_fdsc_sq) > 1e-14) {
          OK = false;

          device::IdxArray<NDim> global_idx;
          layout.putSpatialLocationFromMemoryIndexInto(global_idx, indices...);

          std::stringstream ss;
          ss << "ForwDiff sq mismatch at global (";
          for (size_t i = 0; i < NDim; ++i) {
            ss << global_idx[i];
            if (i < NDim - 1) ss << ", ";
          }
          ss << "): expect = " << expect << ", fdSC_sq = " << val_fdsc_sq << "\n";
          sayMPI << ss.str();
        }
      });
      tdd.verify(OK);
    }

    // Test 3: ForwDiff of product of coordinates (if NDim >= 2)
    if constexpr (NDim >= 2) {
      Field<double, NDim> sc1("SC1_prod", toolBox);
      sc1 = x(1_c);
      sc1.updateGhosts();

      Field<double, NDim> sc2("SC2_prod", toolBox);
      sc2 = x(2_c);
      sc2.updateGhosts();

      Field<double, NDim> sc_prod("SC_prod_field", toolBox);
      sc_prod = sc1 * sc2;
      sc_prod.updateGhosts();

      // Test forward diff in direction 1 (index 0)
      Field<double, NDim> fd_prod_dir1("fd_prod_dir1", toolBox);
      fd_prod_dir1 = forwDiff(sc1 * sc2, Tag<1>{});

      bool OK = true;
      auto sc_prod_view = sc_prod.getFullNDHostView();
      auto fd_prod_dir1_view = fd_prod_dir1.getLocalNDHostView();

      constexpr size_t d = 0;

      NDLoop<NDim>(fd_prod_dir1_view, [&](const auto &...indices) {
        device::IdxArray<NDim> idx_base = {(indices + nGhost)...};
        device::IdxArray<NDim> idx_next = idx_base;
        idx_next[d] += 1;

        const double val_current = device::apply([&](auto... i) { return sc_prod_view(i...); }, idx_base);
        const double val_next = device::apply([&](auto... i) { return sc_prod_view(i...); }, idx_next);
        const double val_fd = fd_prod_dir1_view(indices...);

        const double dx = 1.0;
        const double expect = (val_next - val_current) / dx;

        if (std::abs(expect - val_fd) > 1e-14) {
          OK = false;

          device::IdxArray<NDim> global_idx;
          layout.putSpatialLocationFromMemoryIndexInto(global_idx, indices...);

          std::stringstream ss;
          ss << "ForwDiff prod dir1 mismatch at global (";
          for (size_t i = 0; i < NDim; ++i) {
            ss << global_idx[i];
            if (i < NDim - 1) ss << ", ";
          }
          ss << "): expect = " << expect << ", fd = " << val_fd << "\n";
          sayMPI << ss.str();
        }
      });
      tdd.verify(OK);
    }
  }

  /** @brief BC-aware ForwDiff test: antiperiodic in dim 0, constant field.
   *
   * Sets f = constant c on a lattice with Antiperiodic BC in dim 0. The forward difference in
   * dir 1 (= dim 0) should evaluate to 0 everywhere except at the cell adjacent to the high
   * boundary in dim 0, where the wrap value is -c instead of c, giving (-c - c)/dx = -2c/dx.
   * This isolates the antiperiodic ghost-fill behavior — periodic BC would give 0 everywhere.
   */
  template <size_t NDim> struct ForwDiffBCTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> inline void ForwDiffBCTester<NDim>::Test(TDDAssertion &tdd)
  {
    const device::Idx nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    BCSpec<NDim> spec = allPeriodic<NDim>();
    spec[0] = BCType::Antiperiodic;

    const double c = 1.0;
    const double dx = 1.0; // LatticeParameters<double>() default

    Field<double, NDim> f("f_const", toolBox, LatticeParameters<double>(), spec);
    f = c;
    f.updateGhosts();

    Field<double, NDim> df("df", toolBox);
    df = forwDiff(f, Tag<1>{}); // dir 1 → dim 0

    auto layout = toolBox->mLayouts.getConfigSpaceLayout();
    const auto &localSizes = layout.getLocalSizes();
    const auto &localStarts = layout.getLocalStarts();

    auto view = df.getLocalNDHostView();
    bool ok = true;
    NDLoop<NDim>(view, [&](const auto &...indices) {
      device::IdxArray<NDim> idx{indices...};
      // The LAST owned cell in dim 0 (idx[0] == localSize_0 - 1) on the high-boundary rank is
      // where forwDiff reaches the antiperiodic ghost.
      const bool lastInDim0 = (idx[0] == localSizes[0] - 1);
      const bool highBoundaryRank = (localStarts[0] + localSizes[0] == nGrid);
      const double expected = (lastInDim0 && highBoundaryRank) ? (-2.0 * c / dx) : 0.0;
      const double got = view(indices...);
      if (std::abs(got - expected) > 1e-12) {
        ok = false;
        std::stringstream ss;
        ss << "ForwDiffBC mismatch at idx=" << idx << " got=" << got << " expected=" << expected
           << " (lastInDim0=" << lastInDim0 << " highBoundaryRank=" << highBoundaryRank << ")\n";
        sayMPI << ss.str();
      }
    });
    tdd.verify(ok);
  }

} // namespace TempLat

namespace
{
  // 1D MemoryToolBox is rejected at compile time in MPI builds (see memorytoolbox.h static_assert).
#ifndef HAVE_MPI
  TempLat::TDDContainer<TempLat::ForwDiffTester<1>> test1;
#endif
  TempLat::TDDContainer<TempLat::ForwDiffTester<2>> test2;
  TempLat::TDDContainer<TempLat::ForwDiffTester<3>> test3;
  TempLat::TDDContainer<TempLat::ForwDiffTester<4>> test4;
  TempLat::TDDContainer<TempLat::ForwDiffTester<5>> test5;

  TempLat::TDDContainer<TempLat::ForwDiffBCTester<2>> bcTest2;
  TempLat::TDDContainer<TempLat::ForwDiffBCTester<3>> bcTest3;
} // namespace
