#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_HASKIR_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_HASKIR_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2026

#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A concept to detect if the getKIR method is defined.
   *
   * Unit test: ctest -R test-haskir
   **/
  template <class T>
  concept HasKIR = requires(std::decay_t<T> t) { t.getKIR(); };

} // namespace TempLat

#endif
