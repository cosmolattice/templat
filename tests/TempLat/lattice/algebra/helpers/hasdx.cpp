
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/helpers/hasdx.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasDxTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasDxTester::Test(TDDAssertion &tdd)
  {
    struct myTest {
      double getDx() { return 76; }
    };
    struct myTest2 {
      double getdx() { return 777; }
    };
    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(HasDx<myTest> == true);
    tdd.verify(HasDx<myTest2> == false);
    tdd.verify(HasDx<double> == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasDxTester> test;
}
