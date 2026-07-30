
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/util/rangeiteration/tagliteral.h" // for the _c integer-tag literal
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct MultiplyTester {
    static void Test(TDDAssertion &tdd);
  };

  // A minimal "variable" mock: evaluates to a fixed value, and its symbolic derivative w.r.t. itself is 1.
  // d() returns an evaluable unit (MulVar(1.0)) rather than OneType so that N * d() stays an evaluable
  // expression -- N * OneType constant-folds to a bare scalar which would have no .eval().
  class MulVar
  {
  public:
    MulVar(double b) : v(b) {}
    DEVICE_INLINE_FUNCTION auto eval(const double &) const { return v; }
    template <typename U> auto d(const U &) const { return MulVar(1.0); }

  private:
    double v;
  };

  void MultiplyTester::Test(TDDAssertion &tdd)
  {
    class myClass
    {
    public:
      myClass(int b) : a(b) {}

      DEVICE_INLINE_FUNCTION
      auto eval(device::Idx i) const { return a; }

    private:
      int a;
    };

    myClass c(3);
    myClass b(2);

    //    say << mAdd.get(0, NULL) << " " << GetCPPTypeName::get(decltype(mAdd.get(0, NULL))) << "\n";
    //    say << HasEvalMethod<Operators::Add<GetterGetOffset, GetterGetOffset>>::value << "\n";
    tdd.verify(HasEvalMethod<Operators::Multiplication<myClass, myClass>> == true);
    tdd.verify((b * c).eval(0) == 6);

    int e = 3, f = 4;
    tdd.verify(HasEvalMethod<decltype(e * f)> == false);

    // pointless, but shuts up the compiler about unused variables:
    e = e + f;

    // --- Regression test for G4: d/dx (N*x) must be N. ---
    // MultiplicationN::d() currently returns N * mR (the value N*x) instead of N * GetDeriv::get(mR, other),
    // i.e. it drops the chain rule and returns the value rather than the derivative.
    MulVar x(3.0);
    auto threeX = x * 3_c; // MultiplicationN<MulVar, 3>
    tdd.verify(AlmostEqual(threeX.d(x).eval(0), 3.0));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MultiplyTester> test;
}
