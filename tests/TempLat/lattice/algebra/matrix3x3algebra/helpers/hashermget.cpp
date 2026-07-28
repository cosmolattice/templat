
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hashermget.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasHermGetTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasHermGetTester::Test(TDDAssertion &tdd)
  {
    struct MyStruct {
      double HermGet(Tag<0> t) { return 87; };
    };
    struct MyStruct2 {
      double getComp(Tag<0> t) { return 87; };
    };
    struct MyStruct3 {
      double HermGet(Tag<1> t, Tag<1> t2) { return 64; };
    };

    tdd.verify(HasHermGet<MyStruct> == true);
    tdd.verify(HasHermGet<MyStruct2> == false);
    tdd.verify(HasHermGet<MyStruct3> == true);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasHermGetTester> test;
}
