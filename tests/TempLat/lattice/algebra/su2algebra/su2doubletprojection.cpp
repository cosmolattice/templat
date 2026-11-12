/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/su2algebra/su2doubletprojection.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/util/constexpr_for.h"

namespace TempLat
{
  /** @brief doubletProjection(x, c) equals MakeSU2(b, Real(scalar_prod(x, i_sigma(b) * c))). */
  struct SU2DoubletProjectionTester {
    static void Test(TDDAssertion &tdd);
  };

  inline void SU2DoubletProjectionTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    using T = double;
    auto toolBox = MemoryToolBox<NDim>::makeShared(8, 1);
    SpatialCoordinate<NDim> x(toolBox);

    SU2Doublet<T, NDim> phi("phi", toolBox, LatticeParameters<T>());
    SU2Doublet<T, NDim> chi("chi", toolBox, LatticeParameters<T>());
    constexpr_for<0, 4>([&](auto _c) {
      constexpr Tag<(int)decltype(_c)::value> c;
      const double k = decltype(_c)::value + 1.;
      phi(c) = cos(0.3 * k * x(1_c) - 0.2 * x(2_c) + 0.1 * x(3_c));
      chi(c) = sin(0.2 * x(1_c) + 0.4 * k * x(2_c)) + 0.1 * k * x(3_c);
    });

    auto matches = [&](auto &&a, auto &&b) {
      SU2Field<T, NDim> projected("projected", toolBox, LatticeParameters<T>());
      SU2Field<T, NDim> ref("ref", toolBox, LatticeParameters<T>());
      projected = doubletProjection(a, b);
      ref = MakeSU2(i, Real(scalar_prod(a, Constants::i_sigma(i) * b)));
      double dev = 0;
      constexpr_for<0, 4>([&](auto _c) {
        constexpr Tag<(int)decltype(_c)::value> c;
        dev = std::max(dev, (double)average((projected(c) - ref(c)) * (projected(c) - ref(c))));
      });
      return dev < 1e-28;
    };

    tdd.verify(matches(phi, chi));
    tdd.verify(matches(chi, phi));
    tdd.verify(matches(shift<1>(phi), chi));
  }
} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SU2DoubletProjectionTester> test;
}
