
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#include "TempLat/lattice/algebra/helpers/getngrid.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct GetNGridTester {
    static void Test(TDDAssertion &tdd);
  };

  void GetNGridTester::Test(TDDAssertion &tdd)
  {
    struct MyTestOne {
      MyTestOne() : mToolBox(MemoryToolBox<3>::makeShared(32, 1)) {}
      device::memory::host_ptr<MemoryToolBox<3>> getToolBox() const { return mToolBox; }
      device::memory::host_ptr<MemoryToolBox<3>> mToolBox;
    };
    MyTestOne one;

    struct MyTestTwo {
    };
    MyTestTwo two;

    tdd.verify(GetNGrid::get(one) == 32);
    tdd.verify(GetNGrid::get(two) == 1);
    tdd.verify(GetNGrid::getVec(one)[2] == 32);
    tdd.verify(GetNGrid::getVec(two).size() == 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GetNGridTester> test;
}
