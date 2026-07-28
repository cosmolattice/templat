
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/field/helpers/hasastuplecat.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasAsTupleCatTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasAsTupleCatTester::Test(TDDAssertion &tdd)
  {

    struct myTuple {
      myTuple() : tup(std::make_tuple(1, std::make_tuple(34, 0.965))) {}

      auto asTupleCat() { return std::tuple_cat(tup); };
      std::tuple<int, std::tuple<int, double>> tup;
    };

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(HasAsTupleCat<myTuple>::value == true);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasAsTupleCatTester> test;
}
