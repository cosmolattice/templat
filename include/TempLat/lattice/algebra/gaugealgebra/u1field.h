#ifndef COSMOINTERFACE_GAUGEALGEBRA_U1FIELD_H
#define COSMOINTERFACE_GAUGEALGEBRA_U1FIELD_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/complexalgebra/complexfield.h"
#include "TempLat/lattice/algebra/gaugealgebra/u1exponential.h"

namespace TempLat
{
  /** @brief A U(1) group-valued field.
   *
   * A U(1) group element is a complex number of unit modulus, so a U(1) field is simply a ComplexField:
   * the group product is complex multiplication, the inverse is conj() / dagger(), and the exponential
   * map is complexPhase() (see u1exponential.h).
   *
   * Unlike SU2Field there is no unitarize() method. Links stay on the group by construction when they
   * are updated multiplicatively through the exponential map,
   *
   *   U = complexPhase(dt * E) * U;   // U <- e^{i dt E} U, exactly unit modulus
   *
   * rather than additively (U += ... would leave the group). See bench-u1_evolution for a complete
   * time-evolution example.
   *
   * Unit test: ctest -R test-u1field
   *
   * @vocab-summary A single $U(1)$ link $A_i$ for fixed $i$. An alias of ComplexField, because a $U(1)$ group
   * element is a unit-modulus complex number: the group product is complex multiplication and the inverse is
   *   conj().
   * @vocab-signature U1Field<T, NDim> u("u", toolBox);
   * @vocab-tags ComplexField
   **/
  template <typename T, size_t NDim> using U1Field = ComplexField<T, NDim>;
} // namespace TempLat

#endif
