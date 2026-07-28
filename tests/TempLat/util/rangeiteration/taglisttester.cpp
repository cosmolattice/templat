
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/util/rangeiteration/taglisttester.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/rangeiteration/for_in_range.h"

namespace TempLat
{

  // Okay to define something globally here as not included anywhere

  template <int N> int get(Tag<N> t) { return N; }

  struct TagListTester {
    static void Test(TDDAssertion &tdd);
  };

  template <int N> int get() { return N; }

  void TagListTester::Test(TDDAssertion &tdd)
  {

    int c = 1;
    for_in_range<1, 25>([&](auto i) { tdd.verify(get(i) == c++); });
    c = 1;
    for_in_range<1, 25>([&](auto i) {
      tdd.verify(get<3 * i>() == 3 * c);
      c++;
    });
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::TagListTester> test;
}
