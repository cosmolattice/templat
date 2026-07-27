
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/factorize.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct FactorizeTester {
    static void Test(TDDAssertion &tdd);
  };

  void FactorizeTester::Test(TDDAssertion &tdd)
  {

    tdd.verify(DoesNotThrow<FactorizationFailed>([]() {
      Factorize fac(13);
      Factorize fac1(4);
      Factorize fac2(6);
      Factorize fac3(12);
      Factorize fac4(148);
      if (TDDRegister::isSingleUnitTest()) {
        Factorize fac5(1481451232);
        Factorize fac6(1481444);
        Factorize fac7(1482155);
      }
    }));
    //    std::cerr << fac4 << "\n";
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::FactorizeTester> test;
}
