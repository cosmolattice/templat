/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/algebra/complexalgebra/ascomplexfield.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/complexalgebra/real.h"
#include "TempLat/lattice/algebra/complexalgebra/imag.h"
#include "TempLat/lattice/algebra/complexalgebra/asfourier.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfield.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/lattice/field/field.h"

namespace TempLat
{

  struct AsComplexFieldTester {
    static void Test(TDDAssertion &tdd);
  };

  /** @brief A deterministic complex-valued plain getter, so eval values are checkable exactly. */
  struct MockComplexGetter {
    static constexpr size_t NDim = 3;
    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION complex<double> eval(const IDX &...idx) const
    {
      const device::array<size_t, NDim> c{static_cast<size_t>(idx)...};
      const double re = 1.0 * c[0] + 10.0 * c[1] + 100.0 * c[2];
      const double im = 7.0 * c[0] + 3.0 * c[1] + 1.0 * c[2] + 0.5;
      return complex<double>(re, im);
    }
    std::string toString() const { return "mock"; }
  };

  void AsComplexFieldTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    using T = double;

    MockComplexGetter mock;

    // --- concept plumbing (mutual exclusivity of the Real/Imag overloads) ---
    tdd.verify(HasComplexEval<MockComplexGetter>, "a complex-valued plain getter is adaptable");
    tdd.verify(!HasComplexEval<ComplexField<double, NDim>>, "complex-field expressions are excluded (overload #1)");
    tdd.verify(!HasComplexEval<complex<double>>, "complex values are excluded (overload #2)");
    tdd.verify(!HasComplexEval<Field<double, NDim>>, "real-valued getters are excluded");
    tdd.verify(HasComplexEval<FourierView<double, NDim>>, "a Fourier view's eval is complex -> adaptable");

    // --- adapter shape ---
    auto acf = asComplexField(mock);
    tdd.verify(HasComplexFieldGet<decltype(acf)>, "asComplexField speaks the complex-field protocol");
    static_assert(std::is_same_v<decltype(acf.eval(0ul, 0ul, 0ul)), device::array<double, 2>>,
                  "eval must return array<C,2>, not complex<C> -- assignment loops index result[0]/result[1]");

    const complex<double> ref = mock.eval(1ul, 2ul, 3ul);
    tdd.verify(AlmostEqual(acf.eval(1ul, 2ul, 3ul)[0], device::real(ref)), "adapter eval real component");
    tdd.verify(AlmostEqual(acf.eval(1ul, 2ul, 3ul)[1], device::imag(ref)), "adapter eval imag component");

    // --- Real / Imag third overload ---
    tdd.verify(AlmostEqual(Real(mock).eval(1ul, 2ul, 3ul), device::real(ref)), "Real(mock) == real part");
    tdd.verify(AlmostEqual(Imag(mock).eval(1ul, 2ul, 3ul), device::imag(ref)), "Imag(mock) == imag part");
    tdd.verify(AlmostEqual(Real(mock).eval(1ul, 2ul, 3ul), Real(asComplexField(mock)).eval(1ul, 2ul, 3ul)),
               "Real(f) == Real(asComplexField(f))");

    // --- scalar algebra accepts the components ---
    tdd.verify(AlmostEqual((2 + Real(mock)).eval(1ul, 2ul, 3ul), 2 + device::real(ref)), "2 + Real(mock)");
    tdd.verify(AlmostEqual((Real(mock) * Imag(mock)).eval(1ul, 2ul, 3ul), device::real(ref) * device::imag(ref)),
               "Real(mock) * Imag(mock)");

    // --- round trip: asFourier o asComplexField == identity ---
    const complex<double> rt = asFourier(asComplexField(mock)).eval(1ul, 2ul, 3ul);
    tdd.verify(AlmostEqual(device::real(rt), device::real(ref)) && AlmostEqual(device::imag(rt), device::imag(ref)),
               "asFourier(asComplexField(x)) == x");

    // --- end to end with a real RNG in Fourier space (generation-counter is the load-bearing assertion) ---
    const device::Idx nGrid = 16, nGhost = 2;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    RandomGaussianField<double, NDim> rng("ascomplexfield-seed", toolBox);

    // Pipeline smoke: the adapter works as an RHS in a real ComplexField assignment.
    ComplexField<double, NDim> cref("cref", toolBox);
    cref.inFourierSpace() = asComplexField(rng);

    // Generation-counter: Real(rng) must advance the RNG exactly ONCE per assignment -- the same amount as a
    // bare assignment. A Complexify(Real(rng), Imag(rng)) design would put two rng copies in the tree and
    // advance it twice, drawing the real and imaginary parts from different generations.
    Field<double, NDim> bare("bare", toolBox);
    Field<double, NDim> re("re", toolBox);

    rng.reset();
    bare.inFourierSpace() = rng; // bare complex draw from generation 0
    const std::string stateAfterBare = rng.saveState();

    rng.reset();
    re.inFourierSpace() = Real(rng); // real part of the generation-0 draw
    const std::string stateAfterReal = rng.saveState();

    tdd.verify(stateAfterReal == stateAfterBare, "Real(rng) advances the RNG generation exactly once per assignment");

    // And the value is the real part of the SAME draw: check the first (always-valid) local Fourier cell.
    {
      auto bareHost = bare.inFourierSpace().getRawHostView();
      auto reHost = re.inFourierSpace().getRawHostView();
      tdd.verify(AlmostEqual(device::real(reHost(0)), device::real(bareHost(0))),
                 "Real(rng) reproduces the real part of the bare draw");
      tdd.verify(AlmostEqual(device::imag(reHost(0)), 0.0), "Real(rng) has zero imaginary part");
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::AsComplexFieldTester> test;
}
