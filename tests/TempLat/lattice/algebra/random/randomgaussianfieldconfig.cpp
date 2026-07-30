/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/util/ndloop.h"

#include <cmath>

namespace TempLat
{

  struct RandomGaussianFieldConfigTester {
    static void Test(TDDAssertion &tdd);
  };

  void RandomGaussianFieldConfigTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    using T = double;
    const device::Idx nGrid = 32, nGhost = 2;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // --- the originally-failing line: a config-space random field composes with scalar algebra ---
    Field<T, NDim> phi("phi", toolBox);
    phi = 2 + RandomGaussianFieldConfig<T, NDim>("seed", toolBox);
    {
      auto view = phi.getLocalNDHostView();
      bool finite = true;
      NDLoop<NDim>(view, [&](const auto... idx) { finite &= std::isfinite(view(idx...)); });
      tdd.verify(finite, "phi = 2 + RandomGaussianFieldConfig runs in configuration space and is finite");
    }

    // --- statistics of the raw white noise: mean ~ 0, variance ~ 1 (deterministic; fixed seed) ---
    RandomGaussianFieldConfig<T, NDim> noise("stats-seed", toolBox);
    Field<T, NDim> a("a", toolBox);
    a = noise;
    {
      auto view = a.getLocalNDHostView();
      double sum = 0.0, sumsq = 0.0;
      long count = 0;
      NDLoop<NDim>(view, [&](const auto... idx) {
        const double v = view(idx...);
        sum += v;
        sumsq += v * v;
        ++count;
      });
      const double mean = sum / count;
      const double var = sumsq / count - mean * mean;
      tdd.verify(std::abs(mean) < 0.1, "config white-noise mean is close to 0");
      tdd.verify(var > 0.8 && var < 1.2, "config white-noise variance is close to 1");
    }

    // --- successive assignments differ; reset() reproduces the first draw ---
    Field<T, NDim> b("b", toolBox);
    b = noise; // generation advanced once by the previous assignment to a
    {
      auto va = a.getLocalNDHostView();
      auto vb = b.getLocalNDHostView();
      bool anyDifferent = false;
      NDLoop<NDim>(va, [&](const auto... idx) { anyDifferent |= !AlmostEqual(va(idx...), vb(idx...)); });
      tdd.verify(anyDifferent, "successive draws differ");
    }

    noise.reset();
    Field<T, NDim> c("c", toolBox);
    c = noise;
    {
      auto va = a.getLocalNDHostView();
      auto vc = c.getLocalNDHostView();
      bool reproduced = true;
      NDLoop<NDim>(va, [&](const auto... idx) { reproduced &= AlmostEqual(va(idx...), vc(idx...)); });
      tdd.verify(reproduced, "reset() reproduces the first draw");
    }

    // --- same seed reproduces bit-for-bit; different seeds differ ---
    {
      RandomGaussianFieldConfig<T, NDim> same("stats-seed", toolBox);
      RandomGaussianFieldConfig<T, NDim> other("other-seed", toolBox);
      Field<T, NDim> fSame("fSame", toolBox);
      Field<T, NDim> fOther("fOther", toolBox);
      fSame = same;
      fOther = other;

      auto vs = fSame.getLocalNDHostView();
      auto vo = fOther.getLocalNDHostView();
      auto va = a.getLocalNDHostView();
      bool sameAsA = true, differsFromOther = false;
      NDLoop<NDim>(va, [&](const auto... idx) {
        sameAsA &= AlmostEqual(va(idx...), vs(idx...));
        differsFromOther |= !AlmostEqual(va(idx...), vo(idx...));
      });
      tdd.verify(sameAsA, "same seed reproduces the same field");
      tdd.verify(differsFromOther, "different seeds give different fields");
    }

    // --- generation advances exactly once per assignment ---
    {
      RandomGaussianFieldConfig<T, NDim> gen("gen-seed", toolBox);
      const std::string s0 = gen.saveState();
      Field<T, NDim> f("gf", toolBox);
      f = gen;
      const std::string s1 = gen.saveState();
      tdd.verify(s0 != s1, "assignment advances the generation");
      gen.reset();
      tdd.verify(gen.saveState() == s0, "reset() returns to the initial generation");
    }

    // --- saveState / loadState round-trip ---
    {
      RandomGaussianFieldConfig<T, NDim> rng("serialize-seed", toolBox);
      Field<T, NDim> f1("sf1", toolBox);
      f1 = rng;
      const std::string saved = rng.saveState();

      RandomGaussianFieldConfig<T, NDim> rng2("unrelated-seed", toolBox);
      rng2.loadState(saved);
      tdd.verify(saved == rng2.saveState(), "loadState restores the serialized state");
    }

    // --- space guard: a configuration-space random field used in Fourier space must throw. This assignment
    //     is type-compatible (real eval -> complex Fourier storage), so without the guard it would silently
    //     draw white noise into Fourier modes. It is caught by DimensionCountRecorder::confirmSpace, which
    //     only fires on the host now that its throw is gated on DEVICE_REGION rather than DEVICE_KOKKOS. ---
    {
      Field<T, NDim> wrong("wrong", toolBox);
      bool threw = false;
      try {
        wrong.inFourierSpace() = RandomGaussianFieldConfig<T, NDim>("config-in-fourier", toolBox);
      } catch (const DimensionCountRecorder_CoordinateSpaceException &) {
        threw = true;
      }
      tdd.verify(threw, "a configuration-space random field assigned in Fourier space throws");
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::RandomGaussianFieldConfigTester> test;
}
