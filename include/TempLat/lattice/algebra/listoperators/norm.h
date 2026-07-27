#ifndef TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_NORM_H
#define TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_NORM_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

#include "TempLat/lattice/algebra/operators/squareroot.h"
#include "TempLat/lattice/algebra/operators/power.h"
#include "TempLat/lattice/algebra/listoperators/total.h"

namespace TempLat
{
  /** @brief A class which computes the norm of a list.
   *
   * Unit test: ctest -R test-norm
   *
   * @vocab-summary Squared norm: the sum of the squares of the components. Defined for collections, complex
   * fields and doublets.
   **/
  template <typename R>
    requires(IsTempLatGettable<0, R> || IsSTDGettable<0, R>)
  auto norm2(const R &r)
  {
    return total(pow<2>(r));
  }
} // namespace TempLat

#endif
