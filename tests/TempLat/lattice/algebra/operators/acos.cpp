
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2024

#include "TempLat/lattice/algebra/operators/acos.h"
#include "TempLat/parallel/devices/kokkos/kokkos.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/almostequal.h"

#include <cmath>

namespace TempLat
{

  struct ACosTester {
    static void Test(TDDAssertion &tdd);
  };

  void ACosTester::Test(TDDAssertion &tdd)
  {
    class myClass
    {
    public:
      myClass(double b) : a(b) {}

      DEVICE_INLINE_FUNCTION
      auto eval(const int &i) const { return a; }

    private:
      double a;
    };

    myClass a(0.3);
    tdd.verify(AlmostEqual(acos(a).eval(0), std::acos(0.3)));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ACosTester> test;
}
