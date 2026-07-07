
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2020
#include "TempLat/lattice/algebra/operators/tanh.h"
#include "TempLat/lattice/algebra/operators/operators.h" // for the operators used by Tanh::d() (pow, divide, ...)
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct TanhTester {
    static void Test(TDDAssertion &tdd);
  };

  // A minimal "variable" mock: evaluates to a fixed value, and its symbolic derivative w.r.t. itself is 1.
  class TanhVar
  {
  public:
    TanhVar(double b) : v(b) {}
    DEVICE_INLINE_FUNCTION auto eval(const double &) const { return v; }
    template <typename U> auto d(const U &) const { return OneType(); }

  private:
    double v;
  };

  void TanhTester::Test(TDDAssertion &tdd)
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
    say << tanh(a).eval(0) << "\n";
    tdd.verify(AlmostEqual(tanh(a).eval(0), std::tanh(3.)));

    // --- Regression test for G4: the symbolic derivative d/dx tanh(x) must be sech^2(x) = 1/cosh(x)^2. ---
    // Tanh::d() currently returns GetDeriv::get(mR, other) / pow<2>(sinh(*this)) == 1/sinh(tanh(x))^2 -- wrong
    // function (sinh vs cosh) and wrong argument (*this vs mR).
    TanhVar x(3.0);
    tdd.verify(AlmostEqual(tanh(x).d(x).eval(0), 1.0 / (std::cosh(3.0) * std::cosh(3.0))));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::TanhTester> test;
}
