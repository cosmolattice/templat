#ifndef TEMPLAT_LATTICE_FIELD_HELPERS_HASASTUPLECAT_H
#define TEMPLAT_LATTICE_FIELD_HELPERS_HASASTUPLECAT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include <type_traits>
#include <utility>

namespace TempLat
{
  /** @brief A class which detects whether or not an object has the hasTuplCat method.
   *
   *
   * Unit test: ctest -R test-hasastuplecat
   **/
  template <class, class = std::void_t<>> struct HasAsTupleCat : std::false_type {
  };

  // specialization recognizes types that do have a nested ::type member:
  template <class T> struct HasAsTupleCat<T, std::void_t<decltype(std::declval<T>().asTupleCat())>> : std::true_type {
  };
} // namespace TempLat

#endif
