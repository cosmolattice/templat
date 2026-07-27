
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct for_in_range_Tester {
    static void Test(TDDAssertion &tdd);
  };

  void for_in_range_Tester::Test(TDDAssertion &tdd)
  {
    int c = 1;
    for_in_range<1, 4>([&](auto i) { tdd.verify(i == c++); });
    ForLoop(i, 0, 0, say << (int)i);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::for_in_range_Tester> test;
}
