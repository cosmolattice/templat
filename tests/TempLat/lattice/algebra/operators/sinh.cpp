
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#include "TempLat/lattice/algebra/operators/sinh.h"
#include "TempLat/lattice/algebra/operators/operators.h" // for the operators used by Sinh::d() (unary minus, etc.)
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct SinhTester {
    static void Test(TDDAssertion &tdd);
  };

  // A minimal "variable" mock: evaluates to a fixed value, and its symbolic derivative w.r.t. itself is 1.
  class SinhVar
  {
  public:
    SinhVar(double b) : v(b) {}
    DEVICE_INLINE_FUNCTION auto eval(const double &) const { return v; }
    template <typename U> auto d(const U &) const { return OneType(); }

  private:
    double v;
  };

  void SinhTester::Test(TDDAssertion &tdd)
  {

    class myClass
    {
    public:
      myClass(int b) : a(b) {}

      DEVICE_INLINE_FUNCTION
      auto eval(const double &i) const { return a; }

    private:
      double a;
    };

    myClass a(3);
    // myClass b(4);
    say << sinh(a).eval(0) << "\n";
    tdd.verify(AlmostEqual(sinh(a).eval(0), std::sinh(3.)));

    // --- Regression test for G4: the symbolic derivative d/dx sinh(x) must be cosh(x). ---
    // Sinh::d() currently returns (exp(*this) + exp(-*this)) / 2 == cosh(sinh(x)), using the whole node (*this)
    // instead of the operand mR.
    SinhVar x(3.0);
    tdd.verify(AlmostEqual(sinh(x).d(x).eval(0), std::cosh(3.0)));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SinhTester> test;
}
