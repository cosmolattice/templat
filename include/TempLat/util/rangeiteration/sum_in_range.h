#ifndef TEMPLAT_UTIL_RANGEITERATION_SUM_IN_RANGE_H
#define TEMPLAT_UTIL_RANGEITERATION_SUM_IN_RANGE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2026

#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/algebra/operators/add.h"

namespace TempLat
{
  template <int begin, typename F, int... INT>
  constexpr auto sum_in_range_impl(F &&func, std::integer_sequence<int, INT...> iseq)
  {
    return (func(Tag<begin + INT>()) + ... + ZeroType());
  }

  /** @brief A function which applies + to the results of func(Tag<i>()) for i in [begin, end).
   *
   * Unit test: ctest -R test-sum_in_range
   *
   * @param func the function to apply to each integer in the range before summing
   **/
  template <int begin, int end, typename F> constexpr auto sum_in_range(F &&func)
  {
    return sum_in_range_impl<begin>(std::forward<F>(func), std::make_integer_sequence<int, end - begin>());
  }
} // namespace TempLat

// NOTE: the empty-range guard is keyed on the (end < beg) range test, NOT on `requires { expr; }`.
// A `requires { expr; }` guard would place a nested Total's inner `[&](auto j){...}` lambda inside a
// requires-expression, which clang rejects ("non-local lambda expression cannot have a capture-default"),
// breaking nested `Total(i, ..., Total(j, ...))` (B2SU2, theta_ijk, ...). Gating on the range keeps the
// empty-range safety (for end < beg the `expr` branch is never instantiated) without that hazard.
/**
 * @vocab-summary Sums `expr` over a compile-time index range — the expression-level $\sum_i$.
 * @vocab-signature Total(i, beg, end, expr)
 * @vocab-example auto S = Total(mu, 1, NDim, Total(nu, 1, NDim, trace(plaq(Us, mu, nu))));
 **/
#define Total(i, beg, end, expr)                                                                                       \
  TempLat::sum_in_range<beg, end + 1>([&](auto i) {                                                                    \
    if constexpr ((end) < (beg)) {                                                                                     \
      return TempLat::ZeroType();                                                                                      \
    } else {                                                                                                           \
      return expr;                                                                                                     \
    }                                                                                                                  \
  })

#endif
