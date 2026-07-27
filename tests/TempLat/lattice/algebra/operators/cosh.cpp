
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#include "TempLat/lattice/algebra/operators/cosh.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct CoshTester {
    static void Test(TDDAssertion &tdd);
  };

  // A minimal "variable" mock: evaluates to a fixed value, and its symbolic derivative w.r.t. itself is 1.
  // Defined at namespace scope because local classes may not have member templates.
  class CoshVar
  {
  public:
    CoshVar(double b) : v(b) {}
    DEVICE_INLINE_FUNCTION auto eval(const double &) const { return v; }
    template <typename U> auto d(const U &) const { return OneType(); }

  private:
    double v;
  };

  void CoshTester::Test(TDDAssertion &tdd)
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
    say << cosh(a).eval(0) << "\n";
    tdd.verify(AlmostEqual(cosh(a).eval(0), std::cosh(3.)));

    // --- Regression test for G4: the symbolic derivative d/dx cosh(x) must be sinh(x). ---
    // Cosh::d() currently returns GetDeriv::get(mR, other) * sinh(*this) == sinh(cosh(x)), using the whole
    // node (*this) instead of the operand mR.
    CoshVar x(3.0);
    tdd.verify(AlmostEqual(cosh(x).d(x).eval(0), std::sinh(3.0)));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::CoshTester> test;
}
