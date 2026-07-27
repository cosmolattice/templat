#ifndef TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_TOTAL_H
#define TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_TOTAL_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "listadd.h"
#include "TempLat/util/tuple_size.h"

namespace TempLat
{
  template <typename Tuple, typename F, int... INT>
  constexpr auto total_impl(Tuple &&tup, F &&func, std::integer_sequence<int, INT...> iseq)
  {
    return (func(GetComponent::get(tup, Tag<INT>())) + ... + ZeroType());
  }

  /** @brief Sums all the components of a tuple after applying func to them.
   *
   * Unit test: ctest -R test-total
   *
   *  @param tup the tuple to sum over
   *  @param func the function to apply to each component before summing
   *
   * @vocab-summary Sums every component of a tuple-like expression, optionally applying a function to each
   * first. Reduces over the components, not over the lattice.
   * @vocab-signature total(tup)   total(tup, f)
   **/
  template <class Tuple, typename Function> constexpr auto total(Tuple &&tup, Function &&func)
  {
    return total_impl(std::forward<Tuple>(tup), std::forward<Function>(func),
                      std::make_integer_sequence<int, tuple_size<std::decay_t<Tuple>>::value>());
  }

  /** @brief Sums all the components of a tuple.
   *
   * Unit test: ctest -R test-total
   *
   *  @param tup the tuple to sum over
   **/
  template <class Tuple> constexpr auto total(Tuple &&tup)
  {
    return total_impl(
        std::forward<Tuple>(tup), [](auto x) { return x; },
        std::make_integer_sequence<int, tuple_size<std::decay_t<Tuple>>::value>());
  }

} // namespace TempLat

#endif
