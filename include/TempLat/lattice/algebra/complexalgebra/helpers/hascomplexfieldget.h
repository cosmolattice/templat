#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_HASCOMPLEXFIELDGET_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_HASCOMPLEXFIELDGET_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{

  /** @brief A concept which checks if instance is a complex field.
   *
   * Unit test: ctest -R test-hascomplexfieldget
   **/
  template <class T>
  concept HasComplexFieldGet = requires(T t, Tag<0> tag) { t.ComplexFieldGet(tag); };
} // namespace TempLat

#endif
