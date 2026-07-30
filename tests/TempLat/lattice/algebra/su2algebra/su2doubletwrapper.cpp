
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/su2algebra/su2doubletwrapper.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{
  struct SU2DoubletWrapperTester {
    static void Test(TDDAssertion &tdd);
  };

  void SU2DoubletWrapperTester::Test(TDDAssertion &tdd)
  {
    // Just to check if all compiles

    SU2DoubletWrapper<double, double, double, double> w1(0, 0, 0, 0);
    tdd.verify(w1.toString() == "SU2(0,0,0,0)");

    SU2DoubletWrapper<double, double, double, double> w2(1.0, 2.0, 3.0, 4.0);
    tdd.verify(w2.toString() == "SU2(1,2,3,4)");
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SU2DoubletWrapperTester> test;
}
