#ifndef TEMPLAT_UTIL_MAKELISTTAG_H
#define TEMPLAT_UTIL_MAKELISTTAG_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/util/assignabletuple.h"

namespace TempLat
{

  /** @brief A class which
   *
   *
   * Unit test: ctest -R test-make_list_tag
   **/
  template <int Start, typename F, int... I>
  inline auto make_list_tag_impl(F &&f, std::integer_sequence<int, I...> iseq)
  {
    return make_list(f(Tag<I + Start>())...);
  }
  template <int Start, int End, typename F> inline constexpr auto make_list_tag(F &&f)
  {
    if constexpr (End > Start)
      return make_list_tag_impl<Start>(std::forward<F>(f), std::make_integer_sequence<int, End - Start>());
    else
      return std::tuple<>();
  }

  template <int End, typename F> inline constexpr auto make_list_tag(F &&f)
  {
    if constexpr (End > 0)
      return make_list_tag_impl<0>(std::forward<F>(f), std::make_integer_sequence<int, End>());
    else
      return std::tuple<>();
  }

  template <int Start, typename F, int... I>
  inline auto make_vector_tag_impl(F &&f, std::integer_sequence<int, I...> iseq)
  {
    return make_vector(f(Tag<I + Start>())...);
  }
  template <int Start, int End, typename F> inline constexpr auto make_vector_tag(F &&f)
  {
    if constexpr (End > Start)
      return make_vector_tag_impl<Start>(std::forward<F>(f), std::make_integer_sequence<int, End - Start>());
    else
      return std::tuple<>();
  }

  template <int End, typename F> inline constexpr auto make_vector_tag(F &&f)
  {
    if constexpr (End > 0)
      return make_vector_tag_impl<0>(std::forward<F>(f), std::make_integer_sequence<int, End>());
    else
      return std::tuple<>();
  }

/**
 * @vocab-summary Like MakeVector, but produces a plain tuple rather than a vector expression.
 * @vocab-signature MakeArray(i, beg, end, expr)
 **/
#define MakeArray(i, beg, end, expr)                                                                                   \
  make_list_tag<beg, end + 1>([&](auto i) {                                                                            \
    {                                                                                                                  \
      return expr;                                                                                                     \
    };                                                                                                                 \
  })
/**
 * @vocab-summary Builds a vector expression by instantiating `expr` once for each index in $[beg, end]$, with
 * `i` bound to a compile-time tag.
 * @vocab-signature MakeVector(i, beg, end, expr)
 * @vocab-example auto Bs = MakeVector(i, 1, 3, magneticField(As, i));
 **/
#define MakeVector(i, beg, end, expr) make_vector_tag<beg, end + 1>([&](auto i) { return expr; })

} // namespace TempLat

#endif
