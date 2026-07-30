
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/namedtmpfile.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct NamedTmpFileTester {
    static void Test(TDDAssertion &tdd);
  };

  void NamedTmpFileTester::Test(TDDAssertion &tdd)
  {

    NamedTmpFile ntf;

    tdd.verify(ntf.outfile.tellp() == 0);

    ntf << "Hoi!";

    tdd.verify(ntf.outfile.tellp() != 0);
    tdd.verify(ntf.remove() == 0);
    tdd.verify(ntf.remove() != 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::NamedTmpFileTester> test;
}
