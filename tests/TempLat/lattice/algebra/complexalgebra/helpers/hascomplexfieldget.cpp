
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasComplexFieldGetTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasComplexFieldGetTester::Test(TDDAssertion &tdd)
  {
    struct MyStruct {
      double ComplexFieldGet(Tag<0> t) { return 87; };
    };
    struct MyStruct2 {
      double getComp(Tag<0> t) { return 87; };
    };

    tdd.verify(HasComplexFieldGet<MyStruct> == true);
    tdd.verify(HasComplexFieldGet<MyStruct2> == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasComplexFieldGetTester> test;
}
