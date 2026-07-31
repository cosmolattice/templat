#ifndef TEMPLAT_LATTICE_MEASUREMENTS_ACCUMULATORTYPE_H
#define TEMPLAT_LATTICE_MEASUREMENTS_ACCUMULATORTYPE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <type_traits>

#include "TempLat/parallel/device.h"

namespace TempLat
{
  namespace Helpers
  {
    template <typename T> struct AccumulatorTypeHelper {
      using type = std::conditional_t<std::is_integral_v<T> || std::is_floating_point_v<T>, double, T>;
    };

    template <typename T> struct AccumulatorTypeHelper<complex<T>> {
      using type = complex<double>;
    };
  } // namespace Helpers

  /** @brief The type in which a lattice reduction should accumulate a quantity of type T.
   *
   * A lattice sum runs over N^3 terms of the same sign, so accumulating in the field's own
   * precision costs roughly log2(N^3) bits: on a 512^3 lattice a single-precision sum of energy
   * densities loses all of its significance. Every reduction therefore evaluates the expression in
   * the field type and accumulates in double precision, converting back only at the call site.
   * Non-arithmetic types are left alone.
   **/
  template <typename T> using AccumulatorType = typename Helpers::AccumulatorTypeHelper<T>::type;
} // namespace TempLat

#endif
