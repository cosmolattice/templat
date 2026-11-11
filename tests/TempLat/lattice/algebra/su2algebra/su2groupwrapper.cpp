
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#include "TempLat/lattice/algebra/su2algebra/su2groupwrapper.h"
#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"

namespace TempLat
{

  struct SU2GroupWrapperTester {
    static void Test(TDDAssertion &tdd);
  };

  void SU2GroupWrapperTester::Test(TDDAssertion &tdd)
  {
    // Just to check if all compiles
    SU2GroupWrapper<double, double, double> w1(0., 0., 0.);
    tdd.verify(w1.toString() == "SU2Group(0,0,0)");

    SU2GroupWrapper<double, double, double> w2(2.0, 3.0, 4.0);
    tdd.verify(w2.toString() == "SU2Group(2,3,4)");

    // toSU2 normalizes all four components and keeps the sign of c0: toSU2(-2 I) = -I.
    {
      constexpr size_t NDim = 3;
      auto toolBox = MemoryToolBox<NDim>::makeShared(8, 1);
      SU2Field<double, NDim> U("U", toolBox);
      SU2Field<double, NDim> V("V", toolBox);
      U(0_c) = -2.0;
      U(1_c) = 0.0;
      U(2_c) = 0.0;
      U(3_c) = 0.0;
      V = toSU2(U);
      Field<double, NDim> c0("c0", toolBox);
      c0 = toSU2(U).SU2Get(0_c); // per-component path

      const auto v0 = V(0_c).getLocalNDHostView();
      const auto v1 = V(1_c).getLocalNDHostView();
      const auto v2 = V(2_c).getLocalNDHostView();
      const auto v3 = V(3_c).getLocalNDHostView();
      const auto w0 = c0.getLocalNDHostView();
      bool all_close = true;
      NDLoop<NDim>(v0, [&](const auto... idx) {
        all_close &= std::abs(v0(idx...) + 1) < 1e-15 && std::abs(v1(idx...)) < 1e-15 &&
                     std::abs(v2(idx...)) < 1e-15 && std::abs(v3(idx...)) < 1e-15 && std::abs(w0(idx...) + 1) < 1e-15;
      });
      tdd.verify(all_close);
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SU2GroupWrapperTester> test;
}
