#ifndef TEMPLAT_LATTICE_ALGEBRA_AXIONALGEBRA_MAGNETICFIELD4_H
#define TEMPLAT_LATTICE_ALGEBRA_AXIONALGEBRA_MAGNETICFIELD4_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"

namespace TempLat
{
  /** \brief A function to get the average between the 4 links of the magnetic field from the gauge potential.
   * Specialised to 3D.
   *
   * Unit test: make test-magneticfield
   *
   * @vocab-summary Site-centred magnetic field: the average over the four plaquettes around a site in each
   * direction. Specialised to three dimensions.
   * @vocab-signature magneticField4(Bs)   magneticField4(Bs, i)
   **/
  template <typename R> auto magneticField4(R Bs, Tag<1> t)
  {
    return 0.25f * (Bs + shift(Bs, -2_c) + shift(Bs, -3_c) + shift(shift(Bs, -2_c), -3_c));
  }

  template <typename R> auto magneticField4(R &&Bs, Tag<2> t)
  {
    return 0.25f * (Bs + shift(Bs, -1_c) + shift(Bs, -3_c) + shift(shift(Bs, -1_c), -3_c));
  }

  template <typename R> auto magneticField4(R &&Bs, Tag<3> t)
  {
    return 0.25f * (Bs + shift(Bs, -1_c) + shift(Bs, -2_c) + shift(shift(Bs, -2_c), -1_c));
  }

  template <typename R> auto magneticField4(R &&Bs) { return MakeVector(i, 1, 3, magneticField4(Bs, i)); }

} // namespace TempLat

#endif
