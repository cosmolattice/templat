/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/su2algebra/su2tracedeficit.h"
#include "TempLat/lattice/algebra/su2algebra/su2trace.h"
#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/algebra/su2algebra/su2liealgebrafield.h"
#include "TempLat/lattice/algebra/su2algebra/su2expmap.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/almostequal.h"

#include <cmath>

namespace TempLat
{

  struct SU2TraceDeficitTester {
    static void Test(TDDAssertion &tdd);
  };

  void SU2TraceDeficitTester::Test(TDDAssertion &tdd)
  {
    static constexpr size_t NDim = 3;
    auto toolBox = MemoryToolBox<NDim>::makeShared(8, 1);

    // ------------------------------------------------------------------
    // (1) Identity: 2 - tr(1) = 0, and it must be exactly zero, not roundoff.
    // ------------------------------------------------------------------
    {
      SU2Field<double, NDim> U("U_identity", toolBox);
      U(0_c) = 1.0;
      U(1_c) = 0.0;
      U(2_c) = 0.0;
      U(3_c) = 0.0;

      Field<double, NDim> deficit("deficit_identity", toolBox);
      deficit = traceDeficit(U);

      tdd.verify(deficit.getLocalNDHostView()(0, 0, 0) == 0.0);
    }

    // ------------------------------------------------------------------
    // (2) Agreement with the naive form on both branches, in double, for
    //     generic unit-norm elements where neither form is ill conditioned.
    //     c0 = +0.6 exercises the deficit branch, c0 = -0.6 the direct one.
    // ------------------------------------------------------------------
    for (const double c0 : {0.6, -0.6}) {
      // Distribute the remaining weight over the three algebra components.
      const double ca = std::sqrt((1.0 - c0 * c0) / 3.0);

      SU2Field<double, NDim> U("U_generic", toolBox);
      U(0_c) = c0;
      U(1_c) = ca;
      U(2_c) = ca;
      U(3_c) = ca;

      Field<double, NDim> deficit("deficit_generic", toolBox);
      Field<double, NDim> naive("naive_generic", toolBox);
      deficit = traceDeficit(U);
      naive = 2.0 - trace(U);

      tdd.verify(AlmostEqual(deficit.getLocalNDHostView()(0, 0, 0), 2.0 - 2.0 * c0));
      tdd.verify(AlmostEqual(deficit.getLocalNDHostView()(0, 0, 0), naive.getLocalNDHostView()(0, 0, 0)));
    }

    // ------------------------------------------------------------------
    // (3) The point of the class: near the identity in SINGLE precision the
    //     naive form loses most of its digits while the deficit form does not.
    //
    //     U = exp(theta n), |theta| = 1e-3  =>  2 - tr U = 2(1 - cos|theta|) ~ 1e-6.
    //     Float spacing at c0 ~ 1 is 6e-8, so `2 - 2 c0` carries ~10% relative
    //     error; the deficit form is built from sin|theta| ~ 1e-3 and keeps ~1e-7.
    // ------------------------------------------------------------------
    {
      const double theta = 1e-3;
      const double component = theta / std::sqrt(3.0);
      const double exact = 2.0 * (1.0 - std::cos(theta));

      LatticeParameters<float> latPar;
      SU2LieAlgebraField<float, NDim> algebra("algebra_small", toolBox, latPar);
      algebra(1_c) = static_cast<float>(component);
      algebra(2_c) = static_cast<float>(component);
      algebra(3_c) = static_cast<float>(component);

      SU2Field<float, NDim> U("U_small", toolBox);
      U = exp(algebra);

      Field<float, NDim> deficit("deficit_small", toolBox);
      Field<float, NDim> naive("naive_small", toolBox);
      deficit = traceDeficit(U);
      naive = 2.0f - trace(U);

      const double deficitError = std::abs(deficit.getLocalNDHostView()(0, 0, 0) / exact - 1.0);
      const double naiveError = std::abs(naive.getLocalNDHostView()(0, 0, 0) / exact - 1.0);

      // The deficit form stays at float-level relative accuracy ...
      tdd.verify(deficitError < 1e-5);
      // ... while the naive one is off by a large fraction of the signal. Asserting a LOWER bound on the
      // naive error keeps this a regression test: if `2 - trace` ever becomes accurate here the premise of
      // the class changed and the test should be revisited.
      tdd.verify(naiveError > 1e-3);
    }

    // ------------------------------------------------------------------
    // (4) Robustness to links that have drifted off the group. Scaling a
    //     unit-norm element by (1 + delta) shifts the naive form by an
    //     ADDITIVE -4 delta -- enough to turn it negative near the identity --
    //     whereas the deficit form only moves at relative order delta.
    // ------------------------------------------------------------------
    {
      const double theta = 1e-3;
      const double delta = 1e-5;
      const double component = theta / std::sqrt(3.0);
      const double exact = 2.0 * (1.0 - std::cos(theta));

      SU2Field<double, NDim> U("U_drifted", toolBox);
      U(0_c) = (1.0 + delta) * std::cos(theta);
      U(1_c) = (1.0 + delta) * std::sin(theta) * component / theta;
      U(2_c) = (1.0 + delta) * std::sin(theta) * component / theta;
      U(3_c) = (1.0 + delta) * std::sin(theta) * component / theta;

      Field<double, NDim> deficit("deficit_drifted", toolBox);
      Field<double, NDim> naive("naive_drifted", toolBox);
      deficit = traceDeficit(U);
      naive = 2.0 - trace(U);

      // Deficit form: still within a few times delta of the exact answer, and positive.
      const double deficitValue = deficit.getLocalNDHostView()(0, 0, 0);
      tdd.verify(deficitValue > 0.0);
      tdd.verify(std::abs(deficitValue / exact - 1.0) < 10.0 * delta);

      // Naive form: the -4 delta shift swamps a signal of size 1e-6 and flips the sign.
      tdd.verify(naive.getLocalNDHostView()(0, 0, 0) < 0.0);
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SU2TraceDeficitTester> test;
}
