
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/parameters/parametergetter.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/almostequal.h"

namespace TempLat
{

  struct ParameterGetterTester {
    static void Test(TDDAssertion &tdd);
  };

  void ParameterGetterTester::Test(TDDAssertion &tdd)
  {
    double d = 9.876;

    ParameterGetter<double> pgd(d, "");

    ParameterGetter<int> pgi(d, "xr");

    std::string str = "9.876";

    ParameterGetter<std::string> pgs(str, "str");

    std::string strbis = pgs;

    say << pgs;

    /* Default is to fail: to remind yourself to implement something here. */

    tdd.verify(AlmostEqual(d, pgd()));

    tdd.verify((int)d == pgi());

    tdd.verify(str == pgs());
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ParameterGetterTester> test;
}
