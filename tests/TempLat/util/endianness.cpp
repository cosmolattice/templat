
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/endianness.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct EndiannessTester {
    static void Test(TDDAssertion &tdd);
  };

  void EndiannessTester::Test(TDDAssertion &tdd)
  {

    Endianness endian;
    tdd.verify(endian.isLittle());
    tdd.verify(!endian.isBig());
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::EndiannessTester> test;
}
