#ifndef TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_ZEROTYPE_H
#define TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_ZEROTYPE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"

#include "TempLat/parallel/device.h"

// Forward-declare Tag to break circular dependency with tag.h
namespace TempLat
{
  template <int N> class Tag;
}

namespace TempLat
{
  /** @brief A class which represents zero. Attempt to simplify derivative expressions.
   *
   * Unit test: ctest -R test-zerotype
   *
   * @vocab-summary Compile-time zero. Multiplying by it collapses a whole expression tree to nothing, and
   * adding it is free — this is how TempLat prunes terms before generating any code.
   **/
  struct ZeroType {
    static std::string toString() { return "(ZeroType)0"; }
    template <typename... IDX>
      requires IsVariadicIndex<IDX...>
    DEVICE_INLINE_FUNCTION static constexpr auto eval(const IDX &...i)
    {
      return 0;
    }
    template <int N> constexpr auto operator()(const Tag<N> t) const { return ZeroType(); }
    inline explicit operator double() const { return value; }
    inline explicit operator float() const { return value; }
    static constexpr int value = 0;
  };
} // namespace TempLat

#endif
