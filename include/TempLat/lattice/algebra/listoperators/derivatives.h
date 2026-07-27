#ifndef TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_DERIVATIVES_H
#define TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_DERIVATIVES_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/util/rangeiteration/make_tuple_tag.h"
#include "TempLat/lattice/algebra/helpers/getderiv.h"
#include "TempLat/util/tuple_size.h"

namespace TempLat
{

  /** @brief A class which aaplies vectorially the automatic derivative. Not used anywhere.
   *
   *
   * Unit test: ctest -R test-derivatives
   **/
  class Derivatives
  {
  public:
    // Put public methods here. These should change very little over time.
    Derivatives() = default;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };

  /**
   * @vocab-summary Symbolic derivative of an expression with respect to each field in a collection, returned as
   * a tuple. Provided for completeness; nothing in TempLat uses it.
   **/
  template <typename T, typename R> auto derivatives(T &&expr, R &others)
  {
    constexpr size_t size = tuple_size<R>::value;
    return make_tuple_tag<size>([&](auto i) {
      auto phi = R::Getter::get(others, i);
      return GetDeriv::get(expr, phi);
    });
  }
} // namespace TempLat

#endif
