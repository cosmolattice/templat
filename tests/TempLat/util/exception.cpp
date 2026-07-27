
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/exception.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ExceptionTester {
    static void Test(TDDAssertion &tdd);
  };

  MakeException(TestException);

  void ExceptionTester::Test(TDDAssertion &tdd)
  {

    tdd.verify(Throws<TestException>([] { throw TestException("Hoi!"); }));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ExceptionTester> test;
}
