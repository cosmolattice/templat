#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISARITHMETIC_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISARITHMETIC_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/helpers/iscomplextype.h"
#include <type_traits>

namespace TempLat
{
  template <typename T>
  concept IsArithmetic = (std::is_arithmetic_v<std::decay_t<T>> || IsComplexType<std::decay_t<T>>);
}

#endif