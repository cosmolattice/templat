
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/util/rangeiteration/make_tuple_tag.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/foreach.h"

namespace TempLat
{

  struct make_tuple_tagTester {
    static void Test(TDDAssertion &tdd);
  };

  void make_tuple_tagTester::Test(TDDAssertion &tdd)
  {
    auto test = make_tuple_tag<10>([](auto i) { return 2 * i; });

    int i = 0;
    bool tmp = true;

    for_each(test, [&](auto x) { tmp = tmp && x == 2 * i++; });

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(tmp);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::make_tuple_tagTester> test;
}
