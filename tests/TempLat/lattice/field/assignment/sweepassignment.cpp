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
