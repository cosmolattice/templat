
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/listoperators/listunaryminus.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/foreach.h"
#include "TempLat/util/almostequal.h"

namespace TempLat
{

  struct ListUnaryMinusTester {
    static void Test(TDDAssertion &tdd);
  };

  void ListUnaryMinusTester::Test(TDDAssertion &tdd)
  {

    auto t1 = std::make_tuple(-1., -2., -3., -4.);

    auto t3 = -t1;
    double tmp = 1;

    for_each(t3, [&](auto x) { tdd.verify(x == tmp++); });
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ListUnaryMinusTester> test;
}
