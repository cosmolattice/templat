#ifndef TEMPLAT_UTIL_SHIFTEDINDEXSEQUENCE_H
#define TEMPLAT_UTIL_SHIFTEDINDEXSEQUENCE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include <cstddef>
#include <utility>

namespace TempLat
{
  /** @brief
   * Use as shifted_index_sequence<SHIFT>(std::make_index_sequence<N>{})
   * To get an index sequence<SHIFT,SHIFT+1,...,SHIFT+N-1>
   *
   * Unit test: ctest -R test-shiftedindexsequence
   **/
  template <int SHIFT, size_t... I> auto shifted_index_sequence(std::index_sequence<I...>)
  {
    return std::index_sequence<(I + SHIFT)...>{};
  }
} // namespace TempLat

#endif
