
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelesswrapper.h"

#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct SymTracelessWrapperTester {
    static void Test(TDDAssertion &tdd);
  };

  void SymTracelessWrapperTester::Test(TDDAssertion &tdd)
  {
    /* Default is to fail: to remind yourself to implement something here. */
    auto test = ConstructSymTraceless(1., 2., 3., 4., 5., 4.);

    tdd.verify(AlmostEqual(test.SymTracelessGet(0_c), -2));
    tdd.verify(AlmostEqual(test.SymTracelessGet(1_c), 2));
    tdd.verify(AlmostEqual(test.SymTracelessGet(2_c), 3));
    tdd.verify(AlmostEqual(test.SymTracelessGet(3_c), 1));
    tdd.verify(AlmostEqual(test.SymTracelessGet(4_c), 5));
    tdd.verify(AlmostEqual(test.SymTracelessGet(5_c), 1));

    tdd.verify(AlmostEqual(test.SymTracelessGet(1_c,1_c), -2));
    tdd.verify(AlmostEqual(test.SymTracelessGet(1_c,2_c), 2));
    tdd.verify(AlmostEqual(test.SymTracelessGet(1_c,3_c), 3));
    tdd.verify(AlmostEqual(test.SymTracelessGet(2_c,1_c), 2));
    tdd.verify(AlmostEqual(test.SymTracelessGet(2_c,2_c), 1));
    tdd.verify(AlmostEqual(test.SymTracelessGet(2_c,3_c), 5));
    tdd.verify(AlmostEqual(test.SymTracelessGet(3_c,1_c), 3));
    tdd.verify(AlmostEqual(test.SymTracelessGet(3_c,2_c), 5));
    tdd.verify(AlmostEqual(test.SymTracelessGet(3_c,3_c), 1));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymTracelessWrapperTester> test;
}
