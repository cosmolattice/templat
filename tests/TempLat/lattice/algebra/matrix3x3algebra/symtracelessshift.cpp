
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026, based on test-su2shift by Adrien Florio, 2019.

#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/helpers/getvectorcomponent.h"

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"

namespace TempLat
{

  struct SymTracelessShiftTester {
    static void Test(TDDAssertion &tdd);
  };

  void SymTracelessShiftTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 5;
    device::Idx nGrid = 8, nGhost = 2;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->setVerbose();

    SpatialCoordinate<NDim> x(toolBox);

    SymTracelessField<double, NDim> st("st", toolBox, LatticeParameters<double>());
    st(0_c) = x(1_c);
    st(1_c) = x(2_c);
    st(2_c) = x(3_c);
    st(3_c) = x(4_c);
    st(4_c) = x(5_c);

    // Test SymTracelessShifterByOne: shift the symmetric-traceless field along each spatial direction
    // and verify that each component is shifted correctly. We use 5 dimensions for simplicity.
    {
      constexpr_for<0, NDim>([&](auto _dir) {
        constexpr size_t dir = decltype(_dir)::value;

        // For each component (1_c, 2_c, 3_c, 4_c, 5_c):
        //   If the component index matches the shift direction, the difference should be 1
        //   (or -(nGrid-1) at the periodic boundary).
        //   Otherwise, the difference should be 0.
        constexpr_for<0, 5>([&](auto _comp) {
          constexpr size_t comp = decltype(_comp)::value;

          Field<double, NDim> result("result", toolBox);
          result = shift(st, _dir + Tag<1>()).SymTracelessGet(Tag<comp>()) - st.SymTracelessGet(Tag<comp>());
          auto view = result.getLocalNDHostView();

          bool all_right = true;
          NDLoop<NDim>(view, [&](const auto... idx) {
            if (comp == dir) {
              bool local = (view(idx...) == 1) || (view(idx...) == -(nGrid - 1));
              all_right &= local;
              if (!local)
                std::cout << "ST component " << comp + 1 << " direction " << dir << " value at "
                          << device::IdxArray<NDim>{idx...} << " is " << view(idx...) << ", should be 1" << std::endl;
            } else {
              bool local = (view(idx...) == 0);
              all_right &= local;
              if (!local)
                std::cout << "ST component " << comp + 1 << " direction " << dir << " value at "
                          << device::IdxArray<NDim>{idx...} << " is " << view(idx...) << ", should be 0" << std::endl;
            }
          });
          tdd.verify(all_right);
        });
      });
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymTracelessShiftTester> test;
}
