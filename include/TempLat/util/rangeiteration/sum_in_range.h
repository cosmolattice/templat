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

#define Total(i, beg, end, expr)                                                                                       \
  TempLat::sum_in_range<beg, end + 1>([&](auto i) {                                                                    \
    if constexpr (requires { expr; }) {                                                                                \
      return expr;                                                                                                     \
    } else if constexpr (end < beg) {                                                                                  \
      return TempLat::ZeroType();                                                                                      \
    } else {                                                                                                           \
      static_assert(end < beg, "Total: end must be less than beg if expr is not valid");                               \
      return TempLat::ZeroType();                                                                                      \
    }                                                                                                                  \
  })

#endif
