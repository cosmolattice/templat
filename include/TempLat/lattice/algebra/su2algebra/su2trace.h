#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2TRACE_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2TRACE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/constants/zerotype.h"

namespace TempLat
{
  /** @brief A function which takes the trace of su2 matrices.
   *
   * Unit test: ctest -R test-su2trace
   *
   * @vocab-summary Trace of an SU(2) matrix, $\mathrm{tr}\,U = 2c_0$ in the four-real-component representation.
   **/
  template <typename R> auto trace(R &&r) { return 2 * r.SU2Get(0_c); }

  constexpr inline auto trace(ZeroType r) { return ZeroType(); }
} // namespace TempLat

#endif
