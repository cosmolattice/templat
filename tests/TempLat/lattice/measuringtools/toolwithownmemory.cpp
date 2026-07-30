
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/measuringtools/toolwithownmemory.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ToolWithOwnMemoryTester {
    static void Test(TDDAssertion &tdd);
  };

  void ToolWithOwnMemoryTester::Test(TDDAssertion &tdd)
  {

    auto toolBox = MemoryToolBox<2>::makeShared(16, 1);

    /*typedef double T;

    Field<T> testField("test", toolBox);

    ToolWithOwnMemory<Field<T>> tester(testField);

    auto workSpace = tester.getFieldForMeasurement("");

    tdd.verify( tester.persistentField.get() == nullptr );

    tester.setPersistentMemory();

    auto workSpace2 = tester.getFieldForMeasurement("");

    tdd.verify( tester.persistentField.get() != nullptr && workSpace2.getMemoryManager() ==
    tester.persistentField->getMemoryManager() );*/
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ToolWithOwnMemoryTester> test;
}
