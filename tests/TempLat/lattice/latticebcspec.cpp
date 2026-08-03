
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/lattice/latticebcspec.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{
  struct LatticeBCSpecTester {
    static void Test(TDDAssertion &tdd);
  };

  void LatticeBCSpecTester::Test(TDDAssertion &tdd)
  {
    // Default-constructed LatticeBCSpec is all-periodic — preserves existing behavior
    // for users who don't request a non-default BC.
    LatticeBCSpec<3> def;
    tdd.verify(def.bc == allPeriodic<3>());

    // Explicit construction from an arbitrary BCSpec is byte-equal to the input.
    BCSpec<3> mixed{BCType::Antiperiodic, BCType::Periodic, BCType::Dirichlet};
    LatticeBCSpec<3> lbc(mixed);
    tdd.verify(lbc.bc == mixed);
  }
} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::LatticeBCSpecTester> test;
}
