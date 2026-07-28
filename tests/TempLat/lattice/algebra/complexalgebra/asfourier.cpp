/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/algebra/complexalgebra/asfourier.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/complexalgebra/complexwrapper.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfield.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/field/field.h"

#include <memory>

namespace TempLat
{

  struct AsFourierTester {
    static void Test(TDDAssertion &tdd);
  };

  /** @brief A minimal complex-field expression that records which lifecycle hooks reach it.
   *  Counters are shared_ptr so that copies made while building the expression tree all
   *  increment the same underlying integer. */
  struct LifecycleSpy {
    static constexpr size_t NDim = 3;
    std::shared_ptr<int> pre = std::make_shared<int>(0);
    std::shared_ptr<int> post = std::make_shared<int>(0);
    std::shared_ptr<int> space = std::make_shared<int>(0);

    // Complex-field protocol: real part 3, imaginary part 4.
    double ComplexFieldGet(Tag<0>) const { return 3.0; }
    double ComplexFieldGet(Tag<1>) const { return 4.0; }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...) const
    {
      device::array<double, 2> result;
      result[0] = 3.0;
      result[1] = 4.0;
      return result;
    }

    void preGet() { ++(*pre); }
    void postGet() { ++(*post); }
    template <size_t N> void confirmSpace(const LayoutStruct<N> &, const SpaceStateType &) const { ++(*space); }

    std::string toString() const { return "spy"; }
  };

  void AsFourierTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    using T = double;
    const device::Idx nGrid = 16, nGhost = 2;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // --- lifecycle regression: asFourier must forward pre/postGet and confirmSpace to its subtree ---
    {
      LifecycleSpy spy;
      auto af = asFourier(LifecycleSpy(spy)); // rvalue -> node holds a value copy that shares the counters

      PreGet::apply(af);
      PostGet::apply(af);
      ConfirmSpace::apply(af, toolBox->mLayouts.getFourierSpaceLayout(), SpaceStateType::Fourier);

      tdd.verify(*spy.pre == 1, "asFourier forwards preGet to its subtree");
      tdd.verify(*spy.post == 1, "asFourier forwards postGet to its subtree");
      tdd.verify(*spy.space == 1, "asFourier forwards confirmSpace to its subtree");

      // eval collapses the {re, im} array into a single complex value
      const auto v = af.eval(0ul, 0ul, 0ul);
      tdd.verify(AlmostEqual(device::real(v), 3.0), "asFourier eval real part");
      tdd.verify(AlmostEqual(device::imag(v), 4.0), "asFourier eval imag part");

      // asFourier is itself a plain (complex-valued) getter, not a complex-field expression
      tdd.verify(!HasComplexFieldGet<decltype(af)>, "asFourier yields a plain getter");
    }

    // --- end to end: asFourier reinterprets a (real-component) complex expression as Fourier-space values ---
    {
      Field<T, NDim> out("out", toolBox);

      // Complexify(re, im) has real components; asFourier packs them into complex Fourier data.
      // Smoke test that the assignment runs through the Field pipeline; check the first local cell.
      // (The raw Fourier view has trailing padding, so we do not scan the whole buffer.)
      out.inFourierSpace() = asFourier(Complexify(3.0, 4.0));

      auto outHost = out.inFourierSpace().getRawHostView();
      tdd.verify(AlmostEqual(device::real(outHost(0)), 3.0) && AlmostEqual(device::imag(outHost(0)), 4.0),
                 "asFourier(Complexify(re,im)) packs the components into Fourier-space complex values");
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::AsFourierTester> test;
}
