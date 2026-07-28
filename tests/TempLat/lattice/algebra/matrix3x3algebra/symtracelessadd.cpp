
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"


namespace TempLat
{

  struct SymTracelessAddTester {
    static void Test(TDDAssertion &tdd);
  };


  void SymTracelessAddTester::Test(TDDAssertion &tdd)
  {

    auto symTracelessStruct = ConstructSymTraceless(1., 2., 3., 4., 5., -5.);
    auto symTracelessStruct2 = ConstructSymTraceless(2., 3., 4., 5., 6., -7.);

    auto addtest = symTracelessStruct + symTracelessStruct2;

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(AlmostEqual(addtest(1_c,1_c), 3));
    tdd.verify(AlmostEqual(addtest(1_c,2_c), 5));
    tdd.verify(AlmostEqual(addtest(1_c,3_c), 7));
    tdd.verify(AlmostEqual(addtest(2_c,1_c), 5));
    tdd.verify(AlmostEqual(addtest(2_c,2_c), 9));
    tdd.verify(AlmostEqual(addtest(2_c,3_c), 11));
    tdd.verify(AlmostEqual(addtest(3_c,1_c), 7));
    tdd.verify(AlmostEqual(addtest(3_c,2_c), 11));
    tdd.verify(AlmostEqual(addtest(3_c,3_c), -12));

  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymTracelessAddTester> test;
}
