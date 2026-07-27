
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymtracelessget.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasSymTracelessGetTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasSymTracelessGetTester::Test(TDDAssertion &tdd)
  {
    struct MyStruct {
      double SymTracelessGet(Tag<0> t) { return 87; };
    };
    struct MyStruct2 {
      double getComp(Tag<0> t) { return 87; };
    };
    struct MyStruct3 {
      double SymTracelessGet(Tag<1> t, Tag<1> t2) { return 64; };
    };

    tdd.verify(HasSymTracelessGet<MyStruct> == true);
    tdd.verify(HasSymTracelessGet<MyStruct2> == false);
    tdd.verify(HasSymTracelessGet<MyStruct3> == true);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasSymTracelessGetTester> test;
}
