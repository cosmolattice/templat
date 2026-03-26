
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hasmatrixget.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasMatrixGetTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasMatrixGetTester::Test(TDDAssertion &tdd)
  {
    struct MyStruct {
      double MatrixGet(Tag<0> t) { return 87; };
    };
    struct MyStruct2 {
      double getComp(Tag<0> t) { return 87; };
    };
    struct MyStruct3 {
      double MatrixGet(Tag<1> t, Tag<1> t2) { return 64; };
    };

    tdd.verify(HasMatrixGet<MyStruct> == true);
    tdd.verify(HasMatrixGet<MyStruct2> == false);
    tdd.verify(HasMatrixGet<MyStruct3> == true);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasMatrixGetTester> test;
}
