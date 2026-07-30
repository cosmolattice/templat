
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2doubletget.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasSU2DoubletGetTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasSU2DoubletGetTester::Test(TDDAssertion &tdd)
  {
    struct MyStruct {
      double SU2DoubletGet(Tag<0> t) const { return 87; };
    };
    struct MyStruct2 {
      double SU2Get(Tag<0> t) const { return 87; };
    };

    tdd.verify(HasSU2DoubletGet<MyStruct> == true);
    tdd.verify(HasSU2DoubletGet<MyStruct2> == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasSU2DoubletGetTester> test;
}
