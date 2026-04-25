
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/parallel/device.h"
#include "TempLat/util/ndloop.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/helpers/getvectorcomponent.h"
#include "TempLat/lattice/algebra/operators/subtract.h"

namespace TempLat
{

  template <size_t NDim> struct ExpressionShifterTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> inline void ExpressionShifterTester<NDim>::Test(TDDAssertion &tdd)
  {
    const ptrdiff_t nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    SpatialCoordinate<NDim> x(toolBox);

    // Create fields for each wavenumber component
    std::vector<Field<double, NDim>> phi_components;
    phi_components.reserve(NDim);
    for (size_t d = 0; d < NDim; ++d) {
      phi_components.emplace_back("phi_" + std::to_string(d), toolBox);
    }

    // Set the fields to the components of x
    constexpr_for<0, NDim>([&](auto d) { phi_components[d] = getVectorComponent(x, d); });
    // ExpressionShifter
    {
    }
    // ExpressionShifterByOne
    {
      constexpr_for<0, NDim>([&](auto _i) {
        constexpr size_t i = decltype(_i)::value;
        constexpr_for<0, NDim>([&](auto _j) {
          constexpr size_t j = decltype(_j)::value;
          Field<double, NDim> result("result", toolBox);
          phi_components[i].getMemoryManager()->updateGhosts();
          result = shift(phi_components[i], _j + Tag<1>()) - phi_components[i];
          result.getMemoryManager()->updateGhosts();
          auto view = result.getLocalNDHostView();
          bool all_right = true;

          auto view_orig = phi_components[i].getLocalNDHostView();

          Field<double, NDim> shiftOnly("shiftOnly", toolBox);
          shiftOnly = shift(phi_components[i], _j + Tag<1>());
          shiftOnly.getMemoryManager()->updateGhosts();
          auto view_shift = shiftOnly.getLocalNDHostView();

          NDLoop<NDim>(view, [&](const auto... idx) {
            if (i == j) {
              bool local = (view(idx...) == 1);
              for (size_t k = 0; k < NDim; ++k)
                local = local || (view(idx...) == -(nGrid - 1));

              all_right &= local;

              if (!local) {
                std::cout << "val at " << device::IdxArray<NDim>{idx...} << " is " << view(idx...) << ", should be "
                          << "1 or " << -(nGrid - 1) << "; The original value at this point is " << view_orig(idx...)
                          << " and the shifted value is " << view_shift(idx...) << std::endl;
              }

            } else {
              bool local = (view(idx...) == 0);
              all_right &= local;
              if (!local)
                std::cout << "val at " << device::IdxArray<NDim>{idx...} << " is " << view(idx...) << ", should be "
                          << "0" << std::endl;
            }
          });
          tdd.verify(all_right);
        });
      });
    }
  }

  /** @brief BC-aware shift test: shift a constant field by 1 across the dim-0 boundary, for
   *  each of Antiperiodic, Dirichlet, Neumann. The shifted-by-+1 expression evaluated at
   *  the LAST owned cell on the high-boundary rank reaches the dim-0 high ghost — whose value
   *  is determined by the BC. Other cells just take their right-neighbor's value (= constant).
   */
  template <size_t NDim> struct ShiftBCTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> inline void ShiftBCTester<NDim>::Test(TDDAssertion &tdd)
  {
    const ptrdiff_t nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    auto layout = toolBox->mLayouts.getConfigSpaceLayout();
    const auto &localSizes  = layout.getLocalSizes();
    const auto &localStarts = layout.getLocalStarts();

    const double c = 3.0;

    auto runOne = [&](BCType bc, double expectedAtBoundaryHigh) {
      BCSpec<NDim> spec = allPeriodic<NDim>();
      spec[0] = bc;
      Field<double, NDim> f("f_const_shift", toolBox, LatticeParameters<double>(), spec);
      f = c;
      f.updateGhosts();

      Field<double, NDim> sf("sf", toolBox);
      sf = shift(f, 0_c + Tag<1>{}); // shift by +1 in dim 0

      auto view = sf.getLocalNDHostView();
      bool ok = true;
      NDLoop<NDim>(view, [&](const auto &...indices) {
        device::IdxArray<NDim> idx{indices...};
        const bool lastInDim0 = (idx[0] == localSizes[0] - 1);
        const bool highB = (localStarts[0] + localSizes[0] == nGrid);
        const double expected = (lastInDim0 && highB) ? expectedAtBoundaryHigh : c;
        const double got = view(indices...);
        if (std::abs(got - expected) > 1e-12) {
          ok = false;
          std::stringstream ss;
          ss << "ShiftBC(bc=" << static_cast<int>(bc) << ") mismatch at idx=" << idx
             << " got=" << got << " expected=" << expected << "\n";
          sayMPI << ss.str();
        }
      });
      tdd.verify(ok);
    };

    runOne(BCType::Antiperiodic, -c);   // shift wraps to high-ghost = -c
    runOne(BCType::Dirichlet,    0.0);  // ghost is zero
    runOne(BCType::Neumann,      c);    // ghost mirrors local last owned (= c)
  }

} // namespace TempLat

namespace
{
#ifndef HAVE_MPI
  TempLat::TDDContainer<TempLat::ExpressionShifterTester<1>> test1;
#endif
  TempLat::TDDContainer<TempLat::ExpressionShifterTester<2>> test2;
  TempLat::TDDContainer<TempLat::ExpressionShifterTester<3>> test3;
  TempLat::TDDContainer<TempLat::ExpressionShifterTester<4>> test4;

  TempLat::TDDContainer<TempLat::ShiftBCTester<2>> bcTest2;
  TempLat::TDDContainer<TempLat::ShiftBCTester<3>> bcTest3;
} // namespace
