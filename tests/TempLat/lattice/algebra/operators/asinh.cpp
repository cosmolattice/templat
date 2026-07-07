
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2020
#include "TempLat/lattice/algebra/operators/asinh.h"
#include "TempLat/lattice/algebra/operators/operators.h" // for the operators used by ASinh::d() (sqrt, divide, ...)
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ASinhTester {
    static void Test(TDDAssertion &tdd);
  };

  // A minimal "variable" mock: evaluates to a fixed value, and its symbolic derivative w.r.t. itself is 1.
  class ASinhVar
  {
  public:
    ASinhVar(double b) : v(b) {}
    DEVICE_INLINE_FUNCTION auto eval(const double &) const { return v; }
    template <typename U> auto d(const U &) const { return OneType(); }

  private:
    double v;
  };

  void ASinhTester::Test(TDDAssertion &tdd)
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
    say << asinh(a).eval(0) << "\n";
    tdd.verify(AlmostEqual(asinh(a).eval(0), std::asinh(3.)));

    // --- Regression test for G4: the symbolic derivative d/dx asinh(x) must be 1/sqrt(1 + x^2). ---
    // ASinh::d() currently returns 1 / sqrt(1 + (*this)*(*this)) * ... == 1/sqrt(1 + asinh(x)^2), using the
    // whole node (*this) instead of the operand mR.
    ASinhVar x(3.0);
    tdd.verify(AlmostEqual(asinh(x).d(x).eval(0), 1.0 / std::sqrt(1.0 + 3.0 * 3.0)));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ASinhTester> test;
}
