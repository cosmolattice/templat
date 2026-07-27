
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/algebra/operators/absolutevalue.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct AbsoluteValueTester {
    static void Test(TDDAssertion &tdd);
  };

  void AbsoluteValueTester::Test(TDDAssertion &tdd)
  {
    struct tmpStruct {
      DEVICE_INLINE_FUNCTION
      device::Idx eval(device::Idx) const { return -1; }
    };
    struct tmpStruct2 {
      DEVICE_INLINE_FUNCTION
      complex<double> eval(device::Idx) const { return complex<double>(1, 1); }
    };
    tdd.verify(abs(tmpStruct()).eval(0) == 1);
    tdd.verify(abs(tmpStruct2()).eval(0) == std::sqrt(2));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::AbsoluteValueTester> test;
}
