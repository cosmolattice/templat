#ifndef COSMOINTERFACE_SU2ALGEBRA_HELPERS_HASSU2GET_H
#define COSMOINTERFACE_SU2ALGEBRA_HELPERS_HASSU2GET_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#include "TempLat/util/rangeiteration/tag.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

namespace TempLat
{
  /** @brief A concept which checks whether the object has a su2get method or not.
   *
   * Unit test: ctest -R test-hassu2get
   **/
  template <typename T>
  concept HasSU2Get = requires(T t) { t.SU2Get(0_c); } || requires(T t) { t.SU2Get(1_c); } ||
                      requires(T t) { t.SU2Get(2_c); } || requires(T t) { t.SU2Get(3_c); };
} // namespace TempLat

#endif
