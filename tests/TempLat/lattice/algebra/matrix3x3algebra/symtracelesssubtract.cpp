
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"


namespace TempLat
{

  struct SymTracelessSubtractTester {
    static void Test(TDDAssertion &tdd);
  };


  void SymTracelessSubtractTester::Test(TDDAssertion &tdd)
  {

    auto symTracelessStruct = ConstructSymTraceless(1., 2., 3., 4., 5., -5.);
    auto symTracelessStruct2 = ConstructSymTraceless(2., 4., 6., 8., 10., -10.);

    auto subtracttest = symTracelessStruct - symTracelessStruct2;

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(AlmostEqual(subtracttest(1_c,1_c), -1));
    tdd.verify(AlmostEqual(subtracttest(1_c,2_c), -2));
    tdd.verify(AlmostEqual(subtracttest(1_c,3_c), -3));
    tdd.verify(AlmostEqual(subtracttest(2_c,1_c), -2));
    tdd.verify(AlmostEqual(subtracttest(2_c,2_c), -4));
    tdd.verify(AlmostEqual(subtracttest(2_c,3_c), -5));
    tdd.verify(AlmostEqual(subtracttest(3_c,1_c), -3));
    tdd.verify(AlmostEqual(subtracttest(3_c,2_c), -5));
    tdd.verify(AlmostEqual(subtracttest(3_c,3_c), 5));

  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymTracelessSubtractTester> test;
}
