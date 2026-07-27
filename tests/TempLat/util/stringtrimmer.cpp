
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/stringtrimmer.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct StringTrimmerTester {
    static void Test(TDDAssertion &tdd);
  };

  void StringTrimmerTester::Test(TDDAssertion &tdd)
  {

    std::string full = "   hello world.   ";
    std::string rtd = "   hello world.";
    std::string ltd = "hello world.   ";
    std::string lrtd = "hello world.";

    tdd.verify(StringTrimmer::rtrimmed(full) == rtd);
    tdd.verify(StringTrimmer::ltrimmed(full) == ltd);
    tdd.verify(StringTrimmer::trimmed(full) == lrtd);

    /* verify things weren't in place: */
    tdd.verify(full != rtd);
    tdd.verify(full != ltd);
    tdd.verify(full != lrtd);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::StringTrimmerTester> test;
}
