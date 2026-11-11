
#ifndef TEMPLAT_LATTICE_ALGEBRA_SHIFT_H
#define TEMPLAT_LATTICE_ALGEBRA_SHIFT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/operators/unaryoperator.h"

#include "TempLat/util/tuple_size.h"
#include "TempLat/util/tuple_tools.h"

#include <string>
#include <type_traits>
#include <utility>

namespace TempLat
{
  /** @brief The I-th component of the shift vector (SHIFTS..., 0, 0, ...). */
  template <size_t I, int... SHIFTS> constexpr int shiftComponent()
  {
    if constexpr (I < sizeof...(SHIFTS)) {
      constexpr int s[] = {SHIFTS...};
      return s[I];
    } else
      return 0;
  }

  /** @brief Evaluates r at idx + (SHIFTS...). The indices are made signed first: they may be unsigned (e.g. in
   *  CUDA builds) and the shifts negative. */
  template <int... SHIFTS, typename R, size_t... I, typename... IDX>
  DEVICE_INLINE_FUNCTION auto evalShiftedImpl(const R &r, std::index_sequence<I...>, const IDX &...idx)
  {
    return DoEval::eval(r, (static_cast<std::make_signed_t<IDX>>(idx) +
                            std::integral_constant<int, shiftComponent<I, SHIFTS...>()>::value)...);
  }

  template <int... SHIFTS, typename R, typename... IDX>
  DEVICE_INLINE_FUNCTION auto evalShifted(const R &r, const IDX &...idx)
  {
    return evalShiftedImpl<SHIFTS...>(r, std::index_sequence_for<IDX...>{}, idx...);
  }

  /** @brief reach offset by the shift vector (SHIFTS...), see GhostReach. */
  template <int... SHIFTS> constexpr GhostReach shiftedReach(GhostReach reach)
  {
    size_t axis = 0;
    ((reach = reach.shifted(axis++, SHIFTS)), ...);
    return reach;
  }

  /** @brief The suffix "_(s1,s2,...)" naming shift<s1,s2,...> in expression strings. */
  template <int... SHIFTS> std::string shiftSuffix()
  {
    std::string res = "_(";
    ((res += std::to_string(SHIFTS) + ","), ...);
    res.back() = ')';
    return res;
  }

  /** @brief A class which implements shifts in coordinates.
   *
   * Unit test: ctest -R test-shift
   **/
  template <typename R, int... SHIFTS> class ExpressionShifter : public UnaryOperator<R>
  {
  public:
    using UnaryOperator<R>::mR;

    ExpressionShifter(const R &pR) : UnaryOperator<R>(pR) {}

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return evalShifted<SHIFTS...>(mR, idx...);
    }

    device::Idx confirmGhostsUpToDate(const GhostReach &reach) const
    {
      return ConfirmGhosts::apply(mR, shiftedReach<SHIFTS...>(reach));
    }

    std::string toString() const { return GetString::get(mR) + shiftSuffix<SHIFTS...>(); }
  };

  template <typename R, int _N> class ExpressionShifterByOne : public UnaryOperator<R>
  {
  public:
    static_assert(_N != 0, "_N cannot be 0.");

    static constexpr int N = _N > 0 ? _N : -_N;
    static constexpr int dir = _N > 0 ? 1 : -1;

    using UnaryOperator<R>::mR;

    ExpressionShifterByOne(const R &pR) : UnaryOperator<R>(pR) {}

    template <typename... IDX>
      requires requires(R r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(mR, idx...);
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

  /**
   * @vocab-summary Reads an expression at a site offset by a constant lattice vector: `shift(phi, 1_c)` is
   * $\phi(x+\hat e_1)$. Every finite-difference operator is built from it, and it is what tells TempLat which
   * ghost cells to exchange.
   * @vocab-signature shift(expr, mu)   shift<i,j,k>(expr)
   * @vocab-example auto p = (Us(mu) * shift(Us(nu), mu)) * dagger(shift(Us(mu), nu)) * dagger(Us(nu));
   **/
  template <int... shifts, class R>
    requires((sizeof...(shifts) > 1) && tuple_size<R>::value == 1)
  auto shift(const R &pR)
  {
    return ExpressionShifter<R, shifts...>(pR);
  }

  template <int N, class R>
    requires(tuple_size<R>::value == 1)
  auto shift(const R &pR)
  {
    return ExpressionShifterByOne<R, N>(pR);
  }

  template <class R, int N>
    requires(tuple_size<R>::value == 1)
  auto shift(const R &pR, Tag<N> t)
  {
    return ExpressionShifterByOne<R, N>(pR);
  }

  template <int N> constexpr OneType shift(OneType) { return OneType(); }
  template <int N> constexpr OneType shift(OneType, Tag<N>) { return OneType(); }

  template <int N> constexpr ZeroType shift(ZeroType) { return ZeroType(); }
  template <int N> constexpr ZeroType shift(ZeroType, Tag<N>) { return ZeroType(); }
} // namespace TempLat

#endif
