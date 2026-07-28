#ifndef TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_ONETYPE_H
#define TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_ONETYPE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A class which represents one. Attempt to simplify derivative expressions.
   *
   * Unit test: ctest -R test-zerotype
   *
   * @vocab-summary Compile-time one. Multiplying by it returns the operand unchanged.
   **/
  struct OneType {
    static std::string toString() { return "(OneType)1"; }
    template <typename... IDX>
      requires IsVariadicIndex<IDX...>
    DEVICE_INLINE_FUNCTION static constexpr auto eval(const IDX &...i)
    {
      return 1;
    }
    inline explicit operator double() const { return value; }
    inline explicit operator float() const { return value; }
    static constexpr int value = 1;
  };
} // namespace TempLat

#endif
