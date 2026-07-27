
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2026
//
// Tests the single-argument complex-field wrappers arg(field) and arg2(field) from
// complexalgebra/arg.h. (Previously this lived in complexalgebra/arg.cpp, whose "test-arg" target
// collided with operators/arg.cpp and which did not compile.)
#include "TempLat/lattice/algebra/complexalgebra/arg.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ComplexFieldArgTester {
    static void Test(TDDAssertion &tdd);
  };

  void ComplexFieldArgTester::Test(TDDAssertion &tdd)
  {
    // A minimal evaluable component and a minimal complex-field-like type exposing real/imag parts.
    struct Comp {
      DEVICE_INLINE_FUNCTION double eval(int) const { return val; }
      double val;
    };

    // ComplexFieldGet(Tag<0>) is the real part, ComplexFieldGet(Tag<1>) is the imaginary part.
    struct ComplexField {
      Comp ComplexFieldGet(Tag<0> t) { return {re}; }
      Comp ComplexFieldGet(Tag<1> t) { return {im}; }
      double re, im;
    };

    // z = 3 + 4i  ->  phase = atan2(4, 3) = 0.927...; both arg (in [0, 2pi)) and arg2 (in ]-pi, pi]) agree.
    ComplexField z{3.0, 4.0};
    tdd.verify(AlmostEqual(arg(z).eval(0), 0.9272952180016122));
    tdd.verify(AlmostEqual(arg2(z).eval(0), 0.9272952180016122));

    // z = 3 - 4i  ->  atan2(-4, 3) = -0.927...; arg wraps it to [0, 2pi), arg2 keeps it in ]-pi, pi].
    ComplexField z2{3.0, -4.0};
    tdd.verify(AlmostEqual(arg(z2).eval(0), 5.355890089177974));
    tdd.verify(AlmostEqual(arg2(z2).eval(0), -0.9272952180016122));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ComplexFieldArgTester> test;
}
