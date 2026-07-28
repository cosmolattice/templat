
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/helpers/istemplatgettable.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct IsTempLatGettableTester {
    static void Test(TDDAssertion &tdd);
  };

  void IsTempLatGettableTester::Test(TDDAssertion &tdd)
  {
    /* Default is to fail: to remind yourself to implement something here. */
    struct myDummyStruct {
      int getComp(Tag<0> t) { return 24; }
    };
    struct myDummyStruct2 {
      int get(int t) { return 24; }
    };
    tdd.verify(IsTempLatGettable<0, std::tuple<int, int>> == false);
    tdd.verify(IsTempLatGettable<0, myDummyStruct> == true);
    tdd.verify(IsTempLatGettable<0, myDummyStruct2> == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::IsTempLatGettableTester> test;
}
