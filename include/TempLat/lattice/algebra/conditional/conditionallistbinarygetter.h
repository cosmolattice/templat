#ifndef TEMPLAT_LATTICE_ALGEBRA_CONDITIONAL_CONDITIONALBINARYGETTERVECTOR_H
#define TEMPLAT_LATTICE_ALGEBRA_CONDITIONAL_CONDITIONALBINARYGETTERVECTOR_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/helpers/istemplatgettable.h"
#include "TempLat/lattice/algebra/helpers/isstdgettable.h"

namespace TempLat
{
  /** @brief A class which is part of the implementation of the list algebra. Conditionally returns some binary
   *operation type.
   *
   * Unit test: ctest -R test-conditionalbinarygetterscalarvector
   **/
  template <typename S, typename T>
  concept ConditionalListBinaryGetter =
      (IsSTDGettable<0, S> || IsSTDGettable<0, T> || IsTempLatGettable<0, S> || IsTempLatGettable<0, T>);
} // namespace TempLat

#endif
