#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_HASSYMGET_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_HASSYMGET_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026

#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{

  /** @brief A concept which checks if instance is a symmetric matrix.
   *
   * Unit test: ctest -R test-hassymget
   **/
  template <class T>
  concept HasSymGet1 = requires(T t, Tag<0> tag) { t.SymGet(tag); };

  template <class T>
  concept HasSymGet2 = requires(T t, Tag<1> tag1, Tag<1> tag2) { t.SymGet(tag1, tag2); };

  template <class T>
  concept HasSymGet = HasSymGet1<T> || HasSymGet2<T>;
} // namespace TempLat

#endif
