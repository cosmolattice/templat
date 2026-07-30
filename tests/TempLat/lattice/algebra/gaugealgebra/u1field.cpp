
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/gaugealgebra/u1field.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"

namespace TempLat
{

  struct U1FieldTester {
    static void Test(TDDAssertion &tdd);
  };

  void U1FieldTester::Test(TDDAssertion &tdd)
  {
    static constexpr size_t NDim = 3;
    using T = double;
    const device::Idx nGrid = 16, nGhost = 1;

    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // An angle theta and a "kick" eps, both varying over the lattice.
    SpatialCoordinate x(toolBox);
    Field<T, NDim> theta("theta", toolBox);
    Field<T, NDim> eps("eps", toolBox);
    theta = 0.25 * getVectorComponent(x, 0_c);
    eps = 0.125 * getVectorComponent(x, 1_c);

    // Assignment from the exponential map, then a multiplicative (group) update.
    // The update aliases U on both sides; this is safe because the right hand side is site-local.
    U1Field<T, NDim> U("U", toolBox);
    U = complexPhase(theta);
    U = complexPhase(eps) * U;

    // The result must be e^{i (theta + eps)}, still of unit modulus.
    ComplexField<T, NDim> expected("expected", toolBox);
    expected(0_c) = cos(theta + eps);
    expected(1_c) = sin(theta + eps);

    auto viewR = U(0_c).getLocalNDHostView();
    auto viewI = U(1_c).getLocalNDHostView();
    auto viewRExp = expected(0_c).getLocalNDHostView();
    auto viewIExp = expected(1_c).getLocalNDHostView();

    bool allCorrect = true;
    NDLoop<NDim>(viewR, [&](const auto &...idx) {
      allCorrect = allCorrect && AlmostEqual(viewR(idx...), viewRExp(idx...));
      allCorrect = allCorrect && AlmostEqual(viewI(idx...), viewIExp(idx...));
      allCorrect = allCorrect && AlmostEqual(viewR(idx...) * viewR(idx...) + viewI(idx...) * viewI(idx...), 1.0);
    });
    tdd.verify(allCorrect);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::U1FieldTester> test;
}
