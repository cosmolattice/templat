/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct PostGetTester {
    static void Test(TDDAssertion &tdd);
  };

  class NoPostGet
  {
  public:
    NoPostGet() = default;
  };

  class WithPostGet
  {
  public:
    void postGet() { called = true; }
    static bool called;
  };
  bool WithPostGet::called = false; // Initialize static member variable

  void PostGetTester::Test(TDDAssertion &tdd)
  {
    NoPostGet noPostGet;
    WithPostGet withPostGet;

    PostGet::apply(noPostGet);   // should compile
    PostGet::apply(withPostGet); // should compile

    tdd.verify(WithPostGet::called == true); // should be true, since we called postGet
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::PostGetTester> test;
}