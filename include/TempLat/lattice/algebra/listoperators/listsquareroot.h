#ifndef TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_LISTSQUAREROOT_H
#define TEMPLAT_LATTICE_ALGEBRA_OPERATORS_LISTOPERATORS_LISTSQUAREROOT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

#include "listpower.h"

namespace TempLat
{

  template <typename R>
    requires(IsSTDGettable<0, R> || IsTempLatGettable<0, R>)
  ListPower<R, HalfType> sqrt(const R &r)
  {
    return ListPower<R, HalfType>(r, HalfType());
  }
} // namespace TempLat

#endif
