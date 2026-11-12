#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELDSHIFT_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELDSHIFT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldunaryoperator.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

namespace TempLat
{
  /** @brief A class which implements spatial shifts for complex algebra.
   *
   *
   * Unit test: ctest -R test-complexfieldshift
   **/
  template <typename R, int... N> class ComplexFieldShifter : public ComplexFieldUnaryOperator<R>
  {
  public:
    // Put public methods here. These should change very little over time.

    using ComplexFieldUnaryOperator<R>::mR;

    ComplexFieldShifter(const R &pR) : ComplexFieldUnaryOperator<R>(pR) {}

    template <int M> auto ComplexFieldGet(Tag<M> t) const
    {
      static_assert(M >= 0 && M <= 1, "ComplexFieldGet: M must be 0 or 1 for ComplexFieldShifter");
      return shift<N...>(mR.ComplexFieldGet(t));
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return evalShifted<N...>(mR, idx...);
    }

    device::Idx confirmGhostsUpToDate(const GhostReach &reach) const
    {
      return ConfirmGhosts::apply(mR, shiftedReach<N...>(reach));
    }

    std::string toString() const { return GetString::get(mR) + shiftSuffix<N...>(); }
  };
  template <typename R, int _N> class ComplexFieldShifterByOne : public ComplexFieldUnaryOperator<R>
  {
    static_assert(_N != 0, "_N cannot be 0.");

    static constexpr int N = _N > 0 ? _N : -_N;
    static constexpr int dir = _N > 0 ? 1 : -1;

  public:
    // Put public methods here. These should change very little over time.

    using ComplexFieldUnaryOperator<R>::mR;

    ComplexFieldShifterByOne(const R &pR) : ComplexFieldUnaryOperator<R>(pR) {}

    template <int M> auto ComplexFieldGet(Tag<M> t) const
    {
      static_assert(M >= 0 && M <= 1, "ComplexFieldGet: M must be 0 or 1 for ComplexFieldShifterByOne");
      return shift<_N>(mR.ComplexFieldGet(t));
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return device::apply([&](const auto &...shifted_idx) { return DoEval::eval(mR, shifted_idx...); },
                           tuple_add_to_nth<N - 1, dir>(device::tie(idx...)));
    }

    device::Idx confirmGhostsUpToDate(const GhostReach &reach) const
    {
      return ConfirmGhosts::apply(mR, reach.shifted(N - 1, dir));
    }

    std::string toString() const { return GetString::get(mR) + "_(->" + std::to_string(_N) + ")"; }
  };

  template <int... shifts, class R>
    requires((sizeof...(shifts) > 1) && HasComplexFieldGet<R>)
  auto shift(const R &pR)
  {
    return ComplexFieldShifter<R, shifts...>(pR);
  }

  template <int N, class R>
    requires HasComplexFieldGet<R>
  auto shift(const R &pR)
  {
    return ComplexFieldShifterByOne<R, N>(pR);
  }

  template <class R, int N>
    requires HasComplexFieldGet<R>
  auto shift(const R &pR, Tag<N> t)
  {
    return ComplexFieldShifterByOne<R, N>(pR);
  }
} // namespace TempLat

#endif
