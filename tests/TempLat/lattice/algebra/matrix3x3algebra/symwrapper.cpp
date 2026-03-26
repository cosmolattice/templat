
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/symwrapper.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct SymmetricTracelessWrapperTester {
    static void Test(TDDAssertion &tdd);
  };

  void SymmetricTracelessWrapperTester::Test(TDDAssertion &tdd)
  {
    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(0_c), 1));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(1_c), 2));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(2_c), 3));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(3_c), 4));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(4_c), 5));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(5_c), 6));

    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(1_c,1_c), 1));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(1_c,2_c), 2));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(1_c,3_c), 3));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(2_c,1_c), 2));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(2_c,2_c), 4));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(2_c,3_c), 5));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(3_c,1_c), 3));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(3_c,2_c), 5));
    tdd.verify(AlmostEqual(ConstructSym(1., 2., 3., 4., 5., 6.).SymGet(3_c,3_c), 6));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymmetricTracelessWrapperTester> test;
}
