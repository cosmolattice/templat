#ifndef TEMPLAT_UTIL_MAKEFLATLIST_H
#define TEMPLAT_UTIL_MAKEFLATLIST_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

#include "TempLat/util/flattentuple.h"
#include "TempLat/util/tuplemaker.h"

namespace TempLat
{

  /** @brief A class which creates a flat tuple from potentially a list of composite objects.
   *
   *
   * Unit test: ctest -R test-makeflatlist
   **/
  template <typename... Args> auto make_flat_list(Args... args)
  {
    return flatten_tuple(std::make_tuple(make_tuple_from(args)...));
  }
} // namespace TempLat

#endif
