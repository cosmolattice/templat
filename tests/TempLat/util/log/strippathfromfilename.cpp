
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/log/strippathfromfilename.h"
#include "TempLat/util/tdd/tdd.h"
#include <string>

namespace TempLat
{

  struct StripPathFromFileNameTester {
    static void Test(TDDAssertion &tdd);
  };

  void StripPathFromFileNameTester::Test(TDDAssertion &tdd)
  {

    tdd.verify(StripPathFromFileName("/path/to/hypothetical/file") == std::string("file"));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::StripPathFromFileNameTester> test;
}
