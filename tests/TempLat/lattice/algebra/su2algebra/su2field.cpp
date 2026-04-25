/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/algebra/su2algebra/su2dagger.h"
#include "TempLat/lattice/algebra/su2algebra/su2multiply.h"
#include "TempLat/lattice/algebra/su2algebra/su2shift.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"

namespace TempLat
{
  template <typename T, size_t NDim> struct SU2FieldTester {
    static void Test(TDDAssertion &tdd);
  };

  template <typename T, size_t NDim> void SU2FieldTester<T, NDim>::Test(TDDAssertion &tdd)
  {
    const device::Idx nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    Field<T, NDim> f0("myField0", toolBox);
    Field<T, NDim> f1("myField1", toolBox);
    Field<T, NDim> f2("myField2", toolBox);
    Field<T, NDim> f3("myField3", toolBox);

    auto res = SU2Field<T, NDim>(f0, f1, f2, f3);

    tdd.verify(res.SU2Get(2_c).toString() == "myField2(x)");

    res(1_c) = 6;
    res(2_c) = 12;
    res(3_c) = 24;

    auto ff3 = res(3_c);
    auto ff3_view = ff3.getLocalNDHostView();

    {
      bool all_true = true;
      NDLoop<NDim>(ff3_view, [&](const auto... idx) { all_true &= (ff3_view(idx...) == 24); });
      tdd.verify(all_true);
    }

    SU2Field<double, NDim> mySU2("allNew", toolBox, LatticeParameters<double>());
    tdd.verify(mySU2(3_c).toString() == "allNew_3(x)");
    tdd.verify(mySU2(0_c).toString() == "allNew_0(x)");

    // Verify c0 is initialized to 1 in named constructor
    {
      auto c0 = mySU2(0_c);
      auto c0_view = c0.getLocalNDHostView();
      bool all_one = true;
      NDLoop<NDim>(c0_view, [&](const auto... idx) { all_one &= (c0_view(idx...) == T(1)); });
      tdd.verify(all_one);
    }

    mySU2 = res;
    const auto &fr3 = mySU2(3_c);
    auto fr3_view = fr3.getLocalNDHostView();

    {
      bool all_true = true;
      NDLoop<NDim>(fr3_view, [&](const auto... idx) {
        all_true &= (fr3_view(idx...) == 24);
        if (!all_true) {
          std::cout << "Mismatch at index: ";
          ((std::cout << idx << " "), ...);
          std::cout << " expected 24 but got " << fr3_view(idx...) << std::endl;
        }
      });
      tdd.verify(all_true);
    }

    // Test unitarize: set c1=0.5, c2=0.5, c3=0.5, then unitarize should give c0 = sqrt(1 - 0.75) = 0.5
    {
      SU2Field<double, NDim> uField("uTest", toolBox, LatticeParameters<double>());
      uField(1_c) = 0.5;
      uField(2_c) = 0.5;
      uField(3_c) = 0.5;
      uField.unitarize();

      auto c0 = uField(0_c);
      auto c0_view = c0.getLocalNDHostView();
      bool all_close = true;
      double expected = sqrt(1.0 - 0.75);
      NDLoop<NDim>(c0_view, [&](const auto... idx) { all_close &= (std::abs(c0_view(idx...) - expected) < 1e-14); });
      tdd.verify(all_close);
    }
  }

  /** @brief SU2 antiperiodic-link cross-boundary test.
   *
   * Constructs an SU2Field with Antiperiodic BC in dim 0. Sets every link to U_const = i*σ_1
   * (c0=0, c1=1, c2=c3=0), which is unitary: U U† = σ_1² = I. Computes V(x) = U(x) * U(x+e_0)†.
   * At the cell on the high-boundary rank with x[0] = local_last (i.e., global x[0] = nGrid-1),
   * shift(U, e_0) reads the antiperiodic ghost — which holds U_const† (c1 sign-flipped via the
   * per-component BC dispatch in SU2Field). Then dagger(U_const†) = U_const, so V = U_const *
   * U_const = (iσ_1)² = -I, i.e. (c0,c1,c2,c3) = (-1,0,0,0). Everywhere else V = U_const *
   * U_const† = I = (1,0,0,0). This is the canonical "dagger-on-wrap" sanity check for the
   * antiperiodic gauge link.
   */
  template <size_t NDim> struct SU2AntiperiodicCrossBoundaryTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void SU2AntiperiodicCrossBoundaryTester<NDim>::Test(TDDAssertion &tdd)
  {
    static_assert(NDim >= 2, "SU2 antiperiodic cross-boundary test needs NDim >= 2.");
    const ptrdiff_t nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    BCSpec<NDim> linkBC = allPeriodic<NDim>();
    linkBC[0] = BCType::Antiperiodic;

    SU2Field<double, NDim> U("U", toolBox, linkBC, LatticeParameters<double>());
    U(0_c) = 0.0;
    U(1_c) = 1.0;
    U(2_c) = 0.0;
    U(3_c) = 0.0;
    U.updateGhosts();

    SU2Field<double, NDim> V("V", toolBox);
    V = U * dagger(shift(U, Tag<1>{}));

    auto layout = toolBox->mLayouts.getConfigSpaceLayout();
    const auto &localSizes  = layout.getLocalSizes();
    const auto &localStarts = layout.getLocalStarts();

    auto v0 = V(0_c).getLocalNDHostView();
    auto v1 = V(1_c).getLocalNDHostView();
    auto v2 = V(2_c).getLocalNDHostView();
    auto v3 = V(3_c).getLocalNDHostView();

    bool ok = true;
    NDLoop<NDim>(v0, [&](const auto &...indices) {
      device::IdxArray<NDim> idx{indices...};
      const bool lastInDim0 = (idx[0] == localSizes[0] - 1);
      const bool highB = (localStarts[0] + localSizes[0] == nGrid);
      const bool atBoundary = lastInDim0 && highB;
      const double e0 = atBoundary ? -1.0 : 1.0;
      const double g0 = v0(indices...), g1 = v1(indices...), g2 = v2(indices...), g3 = v3(indices...);
      if (std::abs(g0 - e0) > 1e-12 || std::abs(g1) > 1e-12 ||
          std::abs(g2)      > 1e-12 || std::abs(g3) > 1e-12) {
        ok = false;
        std::stringstream ss;
        ss << "SU2 antiperiodic cross-boundary mismatch at idx=" << idx
           << " got=(" << g0 << "," << g1 << "," << g2 << "," << g3
           << ") expected=(" << e0 << ",0,0,0)\n";
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
  TempLat::TDDContainer<TempLat::SU2FieldTester<double, 1>> test5;
#endif
  TempLat::TDDContainer<TempLat::SU2FieldTester<double, 2>> test6;
  TempLat::TDDContainer<TempLat::SU2FieldTester<double, 3>> test7;
  TempLat::TDDContainer<TempLat::SU2FieldTester<double, 4>> test8;

  TempLat::TDDContainer<TempLat::SU2AntiperiodicCrossBoundaryTester<2>> bcTest2;
  TempLat::TDDContainer<TempLat::SU2AntiperiodicCrossBoundaryTester<3>> bcTest3;
} // namespace
