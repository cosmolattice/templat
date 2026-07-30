
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/parallel/mpi/session/mpiguard.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct MPIGuardTester {
    static void Test(TDDAssertion &tdd);
  };

  void MPIGuardTester::Test(TDDAssertion &tdd)
  {

    if (MPIGuard::GetInstanceCount() < 1) {
      MPIGuard guard(0, NULL, true);
    } else {
      /* there is an instance of MPIGuard in the calling main, which is a good thing. Then we can test if the multiple
       * instantiation protection works. */
      tdd.verify(Throws<MPIGuardInstantiationException>([]() { MPIGuard guard(0, NULL, true); }));
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MPIGuardTester> test;
}
