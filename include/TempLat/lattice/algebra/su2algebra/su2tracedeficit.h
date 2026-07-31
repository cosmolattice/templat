#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2TRACEDEFICIT_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2TRACEDEFICIT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/operators/unaryoperator.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/parallel/device.h"

#include <type_traits>

namespace TempLat
{

  /** @brief The Wilson trace deficit 2 - tr U of an SU(2)-valued expression, returned as a scalar.
   *
   *  Computed WITHOUT the catastrophic cancellation of the naive `2 - trace(U)`. In the four-real-component
   *  representation U = c_0 + i c_a sigma_a with c_0^2 + |c|^2 = 1, the two algebraically equal forms
   *
   *      2 - tr U = 2 (1 - c_0)              (naive)
   *               = 2 |c|^2 / (1 + c_0)      (deficit form)
   *
   *  are conditioned very differently. Near the identity c_0 = 1 - O(|c|^2), so the naive form subtracts two
   *  numbers agreeing to ~2 log10(1/|c|) digits and keeps mostly roundoff, while the deficit form is built
   *  purely from the small components c_a and carries a relative error of order eps. The roles swap as
   *  c_0 -> -1, where 1 + c_0 cancels instead and the naive form is the well-conditioned one - hence the branch.
   *
   *  This matters because the plaquette trace deficit is essentially always consumed with a 1/dx^4 prefactor
   *  (the SU(2) magnetic energy, cf. CosmoInterface FieldFunctionals::B2SU2), which amplifies whatever error
   *  survives. Two effects then compound in single precision:
   *
   *    - the deficit itself is small near the continuum limit, so the naive form loses most of its digits;
   *    - the stored links slowly leave the group (SU2ExpMap is exactly unitary per step, but the accumulated
   *      product U <- exp(dt X) U is never reprojected), giving a plaquette of norm 1 + delta. The naive form
   *      picks that up as an ADDITIVE -4 delta error, comparable to or larger than the signal, whereas the
   *      deficit form only sees it multiplicatively at relative order delta.
   *
   *  Measured on a CosmoLattice SU(2)+Higgs run at N = 64, dx = 0.08 in float: the naive form drives the
   *  magnetic energy NEGATIVE and leaks 6.5e-2 of the total energy over 4.5 time units, while the deficit
   *  form conserves to 6e-6 - a 1e4 improvement from this expression alone, with no reunitarization.
   *
   *  Note the deficit form is manifestly non-negative for any c_0 > -1, whether or not the operand is exactly
   *  unit-norm. That makes it robust against links that have drifted off the group, but it does not repair
   *  such a drift: it reports the deficit of the projected link while the dynamics still uses the drifted one.
   *
   *  Like SU2Dotter, the operand is evaluated ONCE through the fused eval path rather than as four separate
   *  SU2Get(a) component expressions - for a deep operand such as a plaquette or a clover, building the
   *  components separately blows up the expression object.
   *
   * Unit test: ctest -R test-su2tracedeficit
   **/
  template <typename R> class SU2TraceDeficit : public UnaryOperator<R>
  {
  public:
    using UnaryOperator<R>::mR;

    SU2TraceDeficit(const R &pR) : UnaryOperator<R>(pR) {}

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto U = DoEval::eval(mR, idx...);
      using E = std::decay_t<decltype(U[0])>;
      const E c0 = U[0];
      // Use direct comparison rather than a tolerance (device-compatible, cf. PauliVectorsAlgebra).
      // Either branch is well conditioned on its half: 1 + c0 >= 1 on the first, |2 - 2 c0| >= 2 on the second.
      if (c0 > E(0)) return E(2) * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / (E(1) + c0);
      return E(2) - E(2) * c0;
    }

    virtual std::string operatorString() const override { return "2-tr"; }
  };

  /**
   * @vocab-summary Wilson trace deficit $2 - \mathrm{tr}\,U$ of an SU(2) element, evaluated cancellation-free
   * as $2|c|^2/(1+c_0)$ near the identity. Prefer this over `2 - trace(U)` wherever the result is scaled by a
   * large factor, e.g. the plaquette magnetic energy.
   * @vocab-signature traceDeficit(U)
   **/
  template <class R>
    requires HasSU2Get<R>
  auto traceDeficit(const R &r)
  {
    return SU2TraceDeficit<R>(r);
  }

} // namespace TempLat

#endif
