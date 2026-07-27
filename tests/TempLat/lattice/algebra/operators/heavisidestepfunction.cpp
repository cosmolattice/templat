
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/algebra/operators/heavisidestepfunction.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HeavisideStepFunctionTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace TestScratch
  {
    template <int RETURNVALUE> struct HeavisideTesterTemplate {
      DEVICE_INLINE_FUNCTION
      double get(device::Idx i) const { return RETURNVALUE * std::numeric_limits<double>::epsilon(); }
      DEVICE_INLINE_FUNCTION
      double eval(device::Idx i) const { return RETURNVALUE * std::numeric_limits<double>::epsilon(); }
    };
  } // namespace TestScratch

  void HeavisideStepFunctionTester::Test(TDDAssertion &tdd)
  {
    using namespace TestScratch;
    HeavisideTesterTemplate<1> positive;
    HeavisideTesterTemplate<0> zero;
    HeavisideTesterTemplate<-1> negative;

    tdd.verify(heaviside(zero).eval(0) == 1);

    tdd.verify(heaviside(positive).eval(0) == 1);

    tdd.verify(heaviside(negative).eval(0) == 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HeavisideStepFunctionTester> test;
}
