/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/measuringtools/averager.h"

namespace TempLat
{
  /** @brief Expressions mixing fields or coordinates whose memory layouts differ (ghost depth, grid size) throw
   *  LayoutMismatch: the kernel indexes every leaf with the target's memory indices. Distinct toolboxes with equal
   *  layouts mix freely.
   **/
  struct LayoutCompatibilityTester {
    static void Test(TDDAssertion &tdd);
  };

  inline void LayoutCompatibilityTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    auto ghosts1 = MemoryToolBox<NDim>::makeShared(8, 1);
    auto ghosts1Too = MemoryToolBox<NDim>::makeShared(8, 1);
    auto ghosts2 = MemoryToolBox<NDim>::makeShared(8, 2);
    auto grid16 = MemoryToolBox<NDim>::makeShared(16, 1);
    SpatialCoordinate<NDim> x1(ghosts1), x2(ghosts2);
    auto value = [](const auto &x) { return x(1_c) + 10. * x(2_c) + 100. * x(3_c); };

    Field<double, NDim> f1("f1", ghosts1), f1Too("f1Too", ghosts1Too), f2("f2", ghosts2), f16("f16", grid16);
    f1 = value(x1);
    f2 = value(x2);

    // Different ghost depth, either way, pointwise and with a stencil.
    tdd.verify(Throws<LayoutMismatch>([&]() { f2 = f1; }));
    tdd.verify(Throws<LayoutMismatch>([&]() { f1 = f2; }));
    tdd.verify(Throws<LayoutMismatch>([&]() { f2 = LatLapl(f1); }));
    tdd.verify(Throws<LayoutMismatch>([&]() { f2 = f2 + f1; }));
    // Different grid.
    tdd.verify(Throws<LayoutMismatch>([&]() { f16 = f1; }));
    // Coordinates of another layout.
    tdd.verify(Throws<LayoutMismatch>([&]() { f1 = x2(1_c); }));
    // Measurements.
    tdd.verify(Throws<LayoutMismatch>([&]() { average(f1 + f2); }));
    // Composite fields.
    SU2Field<double, NDim> u1("u1", ghosts1), u2("u2", ghosts2);
    tdd.verify(Throws<LayoutMismatch>([&]() { u2 = u1; }));

    // Equal layouts of distinct toolboxes are compatible, and the result is right.
    f1 = value(x1);
    Field<double, NDim> ref("ref", ghosts1);
    ref = 2. * f1 + LatLapl(f1);
    f1Too = 2. * f1 + LatLapl(f1);
    tdd.verify(average(pow<2>(f1Too - ref)) == 0);
    f1Too = value(x1);
    tdd.verify(average(pow<2>(f1Too - f1)) == 0);
  }
} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::LayoutCompatibilityTester> test;
}
