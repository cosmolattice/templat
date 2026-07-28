
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/util/templatvector.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct TempLatVectorTester {
    static void Test(TDDAssertion &tdd);
  };

  void TempLatVectorTester::Test(TDDAssertion &tdd)
  {

    TempLatVector<int> vec(1, 2, 3);

    tdd.verify(vec.getComp(1_c) == 2);

    std::vector<int> sVec(3, 6);

    vec = sVec;

    tdd.verify(vec[2] == 6);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::TempLatVectorTester> test;
}
