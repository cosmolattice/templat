#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISCOMPLEXTYPE_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISCOMPLEXTYPE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2025

#include "TempLat/parallel/device.h"

namespace TempLat
{
  template <typename T>
  concept IsComplexType = std::is_same_v<std::decay_t<T>, complex<typename std::decay_t<T>::value_type>>;
} // namespace TempLat

#endif
