
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/util/istuplelike.h"
#include "TempLat/util/tdd/tdd.h"
#include <array>

namespace TempLat
{

  struct IsTupleLikeTester {
    static void Test(TDDAssertion &tdd);
  };

  void IsTupleLikeTester::Test(TDDAssertion &tdd)
  {

    struct IAmNotATupleAintI {
    };

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(IsTupleLike<std::array<double, 4>>::value == true);
    tdd.verify(IsTupleLike<IAmNotATupleAintI>::value == false);
    tdd.verify(IsTupleLike<int>::value == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::IsTupleLikeTester> test;
}
