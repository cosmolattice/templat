#ifndef TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_SYMBOLS_H
#define TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_SYMBOLS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2024

#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/algebra/constants/onetype.h"
#include "TempLat/lattice/algebra/operators/unaryminus.h"
#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{
  /** @brief A namespace which contains symbols such as the epsilon tensor.
   *
   *
   * Unit test: ctest -R test-symbols
   **/
  namespace Symbols
  {
    template <int I, int J> constexpr DEVICE_FORCEINLINE_FUNCTION auto delta(Tag<I>, Tag<J>) { return ZeroType(); }

    template <int I> constexpr DEVICE_FORCEINLINE_FUNCTION auto delta(Tag<I>, Tag<I>) { return OneType(); }

    // Helper to check for repeated indices
    template <std::size_t N> constexpr bool has_repeats(const std::array<int, N> &arr)
    {
      for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
          if (arr[i] == arr[j]) return true;
      return false;
    }

    // Helper to compute permutation parity (sign)
    template <std::size_t N> consteval int permutation_sign(const std::array<int, N> &arr)
    {
      int sign = 1;
      for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
          if (arr[i] > arr[j]) sign = -sign;
      return sign;
    }

    // Generalized epsilon for arbitrary dimensions
    template <typename... Tags> constexpr DEVICE_FORCEINLINE_FUNCTION auto epsilon(Tags... tags)
    {
      constexpr std::size_t N = sizeof...(Tags);
      constexpr std::array<int, N> idx = {tag_value(tags)...};
      constexpr int sign = permutation_sign(idx);
      if constexpr (has_repeats(idx))
        return ZeroType();
      else if constexpr (sign == 1)
        return OneType();
      else
        return -OneType();
    }
  }; // namespace Symbols
} // namespace TempLat

#endif
