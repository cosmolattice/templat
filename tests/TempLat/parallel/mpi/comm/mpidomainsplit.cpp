/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include "TempLat/parallel/mpi/comm/mpidomainsplit.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{
  struct MPIDomainSplitTester {
    static void Test(TDDAssertion &tdd);
  };

  void MPIDomainSplitTester::Test(TDDAssertion &tdd)
  {
    tdd.verify(DoesNotThrow<MPIDomainSplitException>([]() {
      MPIDomainSplit dom(28, 4, 3);
      MPIDomainSplit dom1(28, 4, 2);
      MPIDomainSplit dom2(32, 3, 2);
      MPIDomainSplit dom3(28, 2, 1);
    }));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MPIDomainSplitTester> test;
}
