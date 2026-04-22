
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/constants/number.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstruct.h"
#include "TempLat/util/almostequal.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  // Dummy lattice-like type with NDim > 0, for testing operator+=(lattice expr)
  // Modelled after myTmpStruct in the averager test.
  template <size_t _NDim> struct DummyLatticeExpr {
    static constexpr size_t NDim = _NDim;
    double val;
    DummyLatticeExpr(double v, device::Idx ngr = 16) : val(v), mt(MemoryToolBox<NDim>::makeShared(ngr, 1)) {}

    template <typename... IDX> DEVICE_INLINE_FUNCTION double eval(const IDX &...) const { return val; }

    auto getToolBox() const { return mt; }
    void confirmSpace(const LayoutStruct<NDim> &, const SpaceStateType &) const {}
    std::string toString() const { return "DummyLatticeExpr"; }

    device::memory::host_ptr<MemoryToolBox<NDim>> mt;
  };

  struct NumberTester {
    static void Test(TDDAssertion &tdd);
  };

  void NumberTester::Test(TDDAssertion &tdd)
  {
    // GetNDim returns 0
    static_assert(GetNDim::get<Number<double>>() == 0, "Number<double> must have NDim = 0");

    // HasEvalMethod is satisfied
    static_assert(HasEvalMethod<Number<double>>, "Number<double> must satisfy HasEvalMethod");

    // IsScalarType is satisfied
    static_assert(IsScalarType<Number<double>>, "Number<double> must satisfy IsScalarType");

    // eval returns the stored value
    Number<double> n{3.14};
    tdd.verify(n.eval(size_t{0}) == 3.14);

    // operator=(T)
    n = 2.71;
    tdd.verify(n.eval(size_t{0}) == 2.71);

    // operator+=(T) — arithmetic overload
    n = 1.0;
    n += 0.5;
    tdd.verify(n.eval(size_t{0}) == 1.5);

    // operator+=(ZeroType) — no-op
    n = 1.0;
    n += ZeroType();
    tdd.verify(n.eval(size_t{0}) == 1.0);

    // operator+=(0-dim expression) — another Number
    Number<double> m{0.25};
    n = 1.0;
    n += m;
    tdd.verify(n.eval(size_t{0}) == 1.25);

    // DoEval works with Number
    Number<double> q{7.0};
    auto result = DoEval::eval(q, size_t{0});
    tdd.verify(result == 7.0);

    // operator= assigns and returns correct value
    Number<double> r{0.0};
    r = 42.0;
    tdd.verify(r.eval(size_t{0}) == 42.0);
    r = -1.5;
    tdd.verify(r.eval(size_t{0}) == -1.5);
    r = 0.0;
    tdd.verify(r.eval(size_t{0}) == 0.0);

    // Scalar algebra: Number * Number produces a valid expression with GetNDim = 0
    Number<double> a{2.0};
    Number<double> b{3.0};
    auto prod = a * b;
    static_assert(GetNDim::get<decltype(prod)>() == 0, "Number * Number must have NDim = 0");
    tdd.verify(DoEval::eval(prod, size_t{0}) == 6.0);

    // Scalar algebra: Number + Number
    auto sum = a + b;
    static_assert(GetNDim::get<decltype(sum)>() == 0, "Number + Number must have NDim = 0");
    tdd.verify(DoEval::eval(sum, size_t{0}) == 5.0);

    // Scalar algebra: Number - Number
    auto diff = a - b;
    static_assert(GetNDim::get<decltype(diff)>() == 0, "Number - Number must have NDim = 0");
    tdd.verify(DoEval::eval(diff, size_t{0}) == -1.0);

    // Scalar algebra: Number * arithmetic
    auto scaled = a * 4.0;
    static_assert(GetNDim::get<decltype(scaled)>() == 0, "Number * double must have NDim = 0");
    tdd.verify(DoEval::eval(scaled, size_t{0}) == 8.0);

    // Scalar algebra: arithmetic * Number
    auto scaled2 = 4.0 * a;
    static_assert(GetNDim::get<decltype(scaled2)>() == 0, "double * Number must have NDim = 0");
    tdd.verify(DoEval::eval(scaled2, size_t{0}) == 8.0);

    // Scalar algebra: compound expression
    Number<double> c{5.0};
    auto compound = a * b + c;
    static_assert(GetNDim::get<decltype(compound)>() == 0, "compound Number expression must have NDim = 0");
    tdd.verify(DoEval::eval(compound, size_t{0}) == 11.0);

    // operator+= with 0-dim expression (scalar algebra result)
    Number<double> acc{0.0};
    acc += a * b;
    tdd.verify(acc.eval(size_t{0}) == 6.0);

    // operator+=(lattice expr) — calls average(), NDim > 0
    // DummyLatticeExpr returns a constant value at every site, so average == that value
    DummyLatticeExpr<3> latticeExpr(4.5);
    static_assert(GetNDim::get<DummyLatticeExpr<3>>() == 3, "DummyLatticeExpr must have NDim = 3");
    Number<double> s{1.0};
    s += latticeExpr;
    tdd.verify(AlmostEqual(s.eval(size_t{0}), 5.5));

    // operator+=(lattice expr) accumulates correctly
    s += latticeExpr;
    tdd.verify(AlmostEqual(s.eval(size_t{0}), 10.0));

    // --- Number * lattice type produces GetNDim = 3 (acceptance criterion) ---
    DummyLatticeExpr<3> latticeExpr2(2.0);
    Number<double> scalar{3.0};
    auto numTimesLattice = scalar * latticeExpr2;
    static_assert(GetNDim::get<decltype(numTimesLattice)>() == 3, "Number * DummyLatticeExpr<3> must have NDim = 3");
    tdd.verify(AlmostEqual(DoEval::eval(numTimesLattice, size_t{0}, size_t{0}, size_t{0}), 6.0));

    auto latticeTimesNum = latticeExpr2 * scalar;
    static_assert(GetNDim::get<decltype(latticeTimesNum)>() == 3, "DummyLatticeExpr<3> * Number must have NDim = 3");
    tdd.verify(AlmostEqual(DoEval::eval(latticeTimesNum, size_t{0}, size_t{0}, size_t{0}), 6.0));

    // --- eval() with multiple indices (variadic correctness) ---
    Number<double> multi{9.5};
    // 1 index
    tdd.verify(multi.eval(size_t{0}) == 9.5);
    // 2 indices
    tdd.verify(multi.eval(size_t{0}, size_t{1}) == 9.5);
    // 3 indices (typical NDim=3 case)
    tdd.verify(multi.eval(size_t{0}, size_t{1}, size_t{2}) == 9.5);
    // 4 indices
    tdd.verify(multi.eval(size_t{0}, size_t{1}, size_t{2}, size_t{3}) == 9.5);

    // --- Number<float> ---
    static_assert(GetNDim::get<Number<float>>() == 0, "Number<float> must have NDim = 0");
    static_assert(HasEvalMethod<Number<float>>, "Number<float> must satisfy HasEvalMethod");
    static_assert(IsScalarType<Number<float>>, "Number<float> must satisfy IsScalarType");

    Number<float> nf{1.5f};
    tdd.verify(nf.eval(size_t{0}) == 1.5f);

    nf = 2.5f;
    tdd.verify(nf.eval(size_t{0}) == 2.5f);

    nf += 0.5f;
    tdd.verify(nf.eval(size_t{0}) == 3.0f);

    nf += ZeroType();
    tdd.verify(nf.eval(size_t{0}) == 3.0f);

    Number<float> nf2{1.0f};
    nf = 0.0f;
    nf += nf2;
    tdd.verify(nf.eval(size_t{0}) == 1.0f);

    auto fprod = nf * nf2;
    static_assert(GetNDim::get<decltype(fprod)>() == 0, "Number<float> * Number<float> must have NDim = 0");
    tdd.verify(DoEval::eval(fprod, size_t{0}) == 1.0f);

    // --- Interaction with ZeroType in expressions ---
    Number<double> z{5.0};
    [[maybe_unused]] auto zProd = z * ZeroType();
    static_assert(std::is_same_v<decltype(zProd), ZeroType>, "Number * ZeroType must return ZeroType");

    [[maybe_unused]] auto zProd2 = ZeroType() * z;
    static_assert(std::is_same_v<decltype(zProd2), ZeroType>, "ZeroType * Number must return ZeroType");

    // --- Interaction with OneType in expressions ---
    Number<double> o{5.0};
    auto oProd = o * OneType();
    tdd.verify(DoEval::eval(oProd, size_t{0}) == 5.0);

    auto oProd2 = OneType() * o;
    tdd.verify(DoEval::eval(oProd2, size_t{0}) == 5.0);

    // --- Default/zero initialization (aggregate) ---
    Number<double> def{};
    tdd.verify(def.eval(size_t{0}) == 0.0);

    Number<float> deff{};
    tdd.verify(deff.eval(size_t{0}) == 0.0f);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::NumberTester> test;
}
