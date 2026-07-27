
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/hash/keccakhashbareclass.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct KeccakHashBareClassTester {
    static void Test(TDDAssertion &tdd);
  };

  void KeccakHashBareClassTester::Test(TDDAssertion &tdd)
  {

    KeccakHashBareClass khbc;
    khbc.compute("Anything.");
    /* check that we cannot compute twice on the same object. */
    tdd.verify(Throws<KeccakHashBareClassReuseException>([&]() { khbc.compute("Anything."); }));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::KeccakHashBareClassTester> test;
}
