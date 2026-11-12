/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/su2algebra/su2field.h"
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

    // unitarize() projects onto SU(2), c_mu -> c_mu / |c|: it keeps the sign of c0 and stays finite when
    // c1^2 + c2^2 + c3^2 > 1.
    {
      SU2Field<double, NDim> uField("uTest", toolBox, LatticeParameters<double>());
      auto checkUnitarize = [&](std::array<double, 4> in, std::array<double, 4> expected, double tol) {
        uField(0_c) = in[0];
        uField(1_c) = in[1];
        uField(2_c) = in[2];
        uField(3_c) = in[3];
        uField.unitarize();

        const auto v0 = uField(0_c).getLocalNDHostView();
        const auto v1 = uField(1_c).getLocalNDHostView();
        const auto v2 = uField(2_c).getLocalNDHostView();
        const auto v3 = uField(3_c).getLocalNDHostView();
        bool all_close = true;
        NDLoop<NDim>(v0, [&](const auto... idx) {
          const double c[4] = {v0(idx...), v1(idx...), v2(idx...), v3(idx...)};
          for (int m = 0; m < 4; ++m)
            all_close &= std::isfinite(c[m]) && std::abs(c[m] - expected[m]) < tol;
        });
        return all_close;
      };

      tdd.verify(checkUnitarize({-1, 0, 0, 0}, {-1, 0, 0, 0}, 1e-15));
      const double s = 1 + 1e-6;
      tdd.verify(checkUnitarize({-0.5 * s, 0.5 * s, 0.5 * s, 0.5 * s}, {-0.5, 0.5, 0.5, 0.5}, 1e-13));
      tdd.verify(checkUnitarize({0, 1 + 1e-12, 0, 0}, {0, 1, 0, 0}, 1e-14));
    }
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
} // namespace
