
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/filetostring.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/namedtmpfile.h"

namespace TempLat
{

  struct FileToStringTester {
    static void Test(TDDAssertion &tdd);
  };

  void FileToStringTester::Test(TDDAssertion &tdd)
  {

    NamedTmpFile ntf;

    std::string payload = "Hello world!";

    ntf << payload;

    ntf.close();

    FileToString fts(ntf.getName());

    std::string ftsRead = fts;

    tdd.verify(ftsRead == payload);

    tdd.verify(ntf.remove() == 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::FileToStringTester> test;
}
