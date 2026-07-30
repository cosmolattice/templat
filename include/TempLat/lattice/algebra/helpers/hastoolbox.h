#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_HASTOOLBOX_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_HASTOOLBOX_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/memory/memorytoolbox.h"

namespace TempLat
{
  /** @brief A concept which determines at compile time whether an object has a method `getToolBox`.
   * See HasGetMethod.
   * Unit test: ctest -R test-hastoolbox
   **/
  template <class T>
  concept HasToolBox = requires(std::decay_t<T> t) { t.getToolBox(); };
} // namespace TempLat

#endif
