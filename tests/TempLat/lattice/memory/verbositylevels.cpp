
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/memory/verbositylevels.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct VerbosityLevelsTester {
    static void Test(TDDAssertion &tdd);
  };

  void VerbosityLevelsTester::Test(TDDAssertion &tdd)
  {
    VerbosityLevels verbosity;

    verbosity.setAllOff();

    tdd.verify(!verbosity.fieldAssignment);

    verbosity.setAllOn();

    tdd.verify(verbosity.fieldAssignment);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::VerbosityLevelsTester> test;
}
