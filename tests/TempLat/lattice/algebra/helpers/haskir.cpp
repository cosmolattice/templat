
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/helpers/haskir.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasKIRTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasKIRTester::Test(TDDAssertion &tdd)
  {
    struct myTest {
      double getKIR() { return 76; }
    };
    struct myTest2 {
      double getdkkhdiu() { return 777; }
    };

    tdd.verify(HasKIR<myTest> == true);
    tdd.verify(HasKIR<myTest2> == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasKIRTester> test;
}
