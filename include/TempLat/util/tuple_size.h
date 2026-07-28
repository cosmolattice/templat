#ifndef TEMPLAT_UTIL_TUPLE_SIZE_H
#define TEMPLAT_UTIL_TUPLE_SIZE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include <cstddef>
#include <tuple>
#include <type_traits>

namespace TempLat
{
  /** @brief A class which overloads tuple_size for the fcn composite objects.
   *
   *
   * Unit test: ctest -R test-tuple_size
   **/
  template <typename Tuple, class = std::void_t<>> struct tuple_size_helper {
    static const size_t value = 1;
  };

  template <typename Tuple> struct tuple_size_helper<Tuple, std::void_t<decltype(std::tuple_size<Tuple>::value)>> {
    static const size_t value = std::tuple_size<Tuple>::value;
  };

  template <typename Tuple, class = std::void_t<>> struct tuple_size {
    static const size_t value = tuple_size_helper<Tuple>::value;
  };

  template <typename Tuple> struct tuple_size<Tuple, std::void_t<decltype(Tuple::size)>> {
    static const size_t value = Tuple::size;
  };
} // namespace TempLat

#endif
