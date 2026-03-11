#ifndef TEMPLAT_LATTICE_ALGEBRA_GAUGEALGEBRA_CENTEREDFIELDS_H
#define TEMPLAT_LATTICE_ALGEBRA_GAUGEALGEBRA_CENTEREDFIELDS_H
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/operators/shift.h"

namespace TempLat
{
  /** \brief A class which computes centered versions of electric and magnetic fields.
   *
   *
   * Unit test: make test-centeredfields
   **/
  class CenteredFields
  {
  public:
    CenteredFields() = delete;

    template <typename ElectricField, int I> static auto E2(ElectricField E, Tag<I> i)
    {
      return 0.5 * (E(i) + shift(E(i), -i));
    }

    template <typename MagneticField> static auto B4(MagneticField B, Tag<1>)
    {
      return 0.25 * (B(1_c) + shift(B(1_c), -2_c) + shift(B(1_c), -3_c) + shift(shift(B(1_c), -2_c), -3_c));
    }

    template <typename MagneticField> static auto B4(MagneticField B, Tag<2>)
    {
      return 0.25 * (B(2_c) + shift(B(2_c), -1_c) + shift(B(2_c), -3_c) + shift(shift(B(2_c), -1_c), -3_c));
    }

    template <typename MagneticField> static auto B4(MagneticField B, Tag<3>)
    {
      return 0.25 * (B(3_c) + shift(B(3_c), -1_c) + shift(B(3_c), -2_c) + shift(shift(B(3_c), -1_c), -2_c));
    }
  };
} // namespace TempLat

#endif
