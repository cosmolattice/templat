
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/algebra/operators/unaryminus.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/constants/halftype.h"

namespace TempLat
{

  struct UnaryMinusTester {
    static void Test(TDDAssertion &tdd);
  };

  class myClass
  {
  public:
    myClass(int b) : a(b) {}

    template <std::integral... IDX> DEVICE_INLINE_FUNCTION auto eval(const IDX &...i) const { return a; }

  private:
    double a;
  };

  void UnaryMinusTester::Test(TDDAssertion &tdd)
  {
    myClass a(4);
    // myClass b(4);
    tdd.verify(AlmostEqual((-a).eval(0), -4));
    tdd.verify(AlmostEqual((-HalfType()).eval(0, 0, 0), -0.5));

    // Sanity: a double negation evaluates to the original value.
    tdd.verify(AlmostEqual((-(-a)).eval(0), 4));

    // --- Regression test for G4: a triple negation must evaluate to -x. ---
    // The double-minus collapse overload operator-(UnaryMinus<UnaryMinus<T>>&&) returns its argument
    // unchanged, so -(-(-a)) collapses to +a (== 4) instead of -a (== -4).
    tdd.verify(AlmostEqual((-(-(-a))).eval(0), -4));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::UnaryMinusTester> test;
}
