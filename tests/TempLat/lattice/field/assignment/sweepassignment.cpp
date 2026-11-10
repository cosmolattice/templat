/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/subtract.h"
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/util/constexpr_for.h"

namespace TempLat
{
  /** @brief Assignments whose right-hand side reads the target at neighbouring sites give the same result as
   *  assigning into a separate field.
   **/
  struct SweepAssignmentTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace SweepAssignmentDetail
  {
    /** @brief Checks phi = expr(phi) in place against ref = expr(copy of phi), for the given band width
     *  (0: automatic). */
    template <size_t NDim, typename Expr>
    bool matchesOutOfPlace(const device::memory::host_ptr<MemoryToolBox<NDim>> &toolBox, device::Idx bandWidth,
                           Expr &&expr)
    {
      SpatialCoordinate<NDim> x(toolBox);
      Field<double, NDim> phi("phi", toolBox);
      Field<double, NDim> copy("copy", toolBox);
      Field<double, NDim> ref("ref", toolBox);

      // No symmetry, so that a read at the wrong site shows up.
      if constexpr (NDim == 3)
        phi = x(1_c) + 10. * x(2_c) + 100. * x(3_c) + x(1_c) * x(2_c) * x(3_c);
      else
        phi = x(1_c) + 10. * x(2_c) + x(1_c) * x(1_c) * x(2_c);
      copy = phi;

      toolBox->mSweepBuffer.setBandWidth(bandWidth);
      ref = expr(copy);
      phi = expr(phi);
      toolBox->mSweepBuffer.setBandWidth(0);

      return average((phi - ref) * (phi - ref)) == 0;
    }

    /** @brief As matchesOutOfPlace, for a composite field F with components first ... last. */
    template <typename F, int first, int last, typename Expr>
    bool compositeMatchesOutOfPlace(const device::memory::host_ptr<MemoryToolBox<3>> &toolBox,
                                    device::Idx bandWidth, Expr &&expr)
    {
      SpatialCoordinate<3> x(toolBox);
      F f("f", toolBox, LatticeParameters<double>());
      F copy("copy", toolBox, LatticeParameters<double>());
      F ref("ref", toolBox, LatticeParameters<double>());

      constexpr_for<first, last + 1>([&](auto _c) {
        constexpr Tag<(int)decltype(_c)::value> c;
        const double k = decltype(_c)::value + 1.;
        f(c) = k * x(1_c) + 10. * x(2_c) + 100. * x(3_c) + x(1_c) * x(2_c) * x(3_c) / k;
        copy(c) = f(c);
      });

      toolBox->mSweepBuffer.setBandWidth(bandWidth);
      ref = expr(copy);
      f = expr(f);
      toolBox->mSweepBuffer.setBandWidth(0);

      bool matches = true;
      constexpr_for<first, last + 1>([&](auto _c) {
        constexpr Tag<(int)decltype(_c)::value> c;
        matches = matches && average((f(c) - ref(c)) * (f(c) - ref(c))) == 0;
      });
      return matches;
    }
  } // namespace SweepAssignmentDetail

  inline void SweepAssignmentTester::Test(TDDAssertion &tdd)
  {
    using SweepAssignmentDetail::matchesOutOfPlace;
    constexpr size_t NDim = 3;
    auto toolBox = MemoryToolBox<NDim>::makeShared(8, 2);

    auto lapl = [](const auto &f) { return 2. * LatLapl(f); };
    auto lapl2 = [](const auto &f) { return f + 0.1 * LatLapl(LatLapl(f)); };
    auto down = [](const auto &f) { return shift<-2>(f); };
    auto up = [](const auto &f) { return shift<1>(f); };
    auto inPlane = [](const auto &f) { return shift<2>(f) + shift<-3>(f); };
    auto diagonal = [](const auto &f) { return shift<-1, 1, 0>(f) - f; };

    // Automatic band width.
    tdd.verify(matchesOutOfPlace(toolBox, 0, lapl));
    tdd.verify(matchesOutOfPlace(toolBox, 0, lapl2));
    tdd.verify(matchesOutOfPlace(toolBox, 0, down));
    tdd.verify(matchesOutOfPlace(toolBox, 0, up));
    tdd.verify(matchesOutOfPlace(toolBox, 0, inPlane));
    tdd.verify(matchesOutOfPlace(toolBox, 0, diagonal));

    // Band widths: one plane, a divisor of the local extent, a non-divisor (partial last band), and wider than
    // the local extent (a single step).
    for (device::Idx b : {1, 2, 3, 100}) {
      tdd.verify(matchesOutOfPlace(toolBox, b, lapl));
      tdd.verify(matchesOutOfPlace(toolBox, b, lapl2));
      tdd.verify(matchesOutOfPlace(toolBox, b, down));
      tdd.verify(matchesOutOfPlace(toolBox, b, up));
      tdd.verify(matchesOutOfPlace(toolBox, b, diagonal));
    }

    // Two dimensions.
    {
      auto toolBox2 = MemoryToolBox<2>::makeShared(8, 2);
      tdd.verify(matchesOutOfPlace(toolBox2, 0, lapl2));
      tdd.verify(matchesOutOfPlace(toolBox2, 3, down));
    }

    // Composite fields, assigned in a single kernel.
    {
      using SweepAssignmentDetail::compositeMatchesOutOfPlace;
      using C = ComplexField<double, NDim>;
      using U = SU2Field<double, NDim>;
      using D = SU2Doublet<double, NDim>;
      using A = SU2LieAlgebraField<double, NDim>;
      using H = SymTracelessField<double, NDim>;
      auto complexExpr = [](const auto &f) { return shift<-1>(f) * shift<2>(f) + f; };
      auto su2Expr = [](const auto &f) { return shift<-1>(f) * dagger(shift<2>(f)); };
      auto sum = [](const auto &f) { return shift<-1>(f) + shift<1>(f); };
      auto symExpr = [](const auto &f) { return shift<-1>(f) + 2. * shift<2>(f); };
      for (device::Idx b : {0, 3}) {
        tdd.verify(compositeMatchesOutOfPlace<C, 0, 1>(toolBox, b, complexExpr));
        tdd.verify(compositeMatchesOutOfPlace<U, 0, 3>(toolBox, b, su2Expr));
        tdd.verify(compositeMatchesOutOfPlace<D, 0, 3>(toolBox, b, sum));
        tdd.verify(compositeMatchesOutOfPlace<A, 1, 3>(toolBox, b, sum));
        tdd.verify(compositeMatchesOutOfPlace<H, 0, 4>(toolBox, b, symExpr));
      }
    }

    // The buffer is allocated only once an assignment needs it.
    {
      auto fresh = MemoryToolBox<NDim>::makeShared(8, 1);
      Field<double, NDim> phi("phi", fresh);
      Field<double, NDim> chi("chi", fresh);
      phi = 1.;
      chi = LatLapl(phi);
      phi = 2. * phi + chi;
      tdd.verify(!fresh->mSweepBuffer.isAllocated());
      phi = LatLapl(phi);
      tdd.verify(fresh->mSweepBuffer.isAllocated());
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SweepAssignmentTester> test;
}
