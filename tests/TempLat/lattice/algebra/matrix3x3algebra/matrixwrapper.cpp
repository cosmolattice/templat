
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/matrixwrapper.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct SymmetricTracelessWrapperTester {
    static void Test(TDDAssertion &tdd);
  };

  void SymmetricTracelessWrapperTester::Test(TDDAssertion &tdd)
  {

    auto test = ConstructMatrix3x3(1., 2, 3, 4., 5, 6., 7, 8, 9);

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(AlmostEqual(test.MatrixGet(0_c), 1));
    tdd.verify(AlmostEqual(test.MatrixGet(1_c), 2));
    tdd.verify(AlmostEqual(test.MatrixGet(2_c), 3));
    tdd.verify(AlmostEqual(test.MatrixGet(3_c), 4));
    tdd.verify(AlmostEqual(test.MatrixGet(4_c), 5));
    tdd.verify(AlmostEqual(test.MatrixGet(5_c), 6));
    tdd.verify(AlmostEqual(test.MatrixGet(6_c), 7));
    tdd.verify(AlmostEqual(test.MatrixGet(7_c), 8));
    tdd.verify(AlmostEqual(test.MatrixGet(8_c), 9));

    tdd.verify(AlmostEqual(test.MatrixGet(1_c,1_c), 1));
    tdd.verify(AlmostEqual(test.MatrixGet(1_c,2_c), 2));
    tdd.verify(AlmostEqual(test.MatrixGet(1_c,3_c), 3));
    tdd.verify(AlmostEqual(test.MatrixGet(2_c,1_c), 4));
    tdd.verify(AlmostEqual(test.MatrixGet(2_c,2_c), 5));
    tdd.verify(AlmostEqual(test.MatrixGet(2_c,3_c), 6));
    tdd.verify(AlmostEqual(test.MatrixGet(3_c,1_c), 7));
    tdd.verify(AlmostEqual(test.MatrixGet(3_c,2_c), 8));
    tdd.verify(AlmostEqual(test.MatrixGet(3_c,3_c), 9));

  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymmetricTracelessWrapperTester> test;
}
