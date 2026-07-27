
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/containsspace.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ContainsSpaceTester {
    static void Test(TDDAssertion &tdd);
  };

  void ContainsSpaceTester::Test(TDDAssertion &tdd)
  {

    tdd.verify(ContainsSpace::test("This has a space") == true);
    tdd.verify(ContainsSpace::test("Thishasnospace") == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ContainsSpaceTester> test;
}
