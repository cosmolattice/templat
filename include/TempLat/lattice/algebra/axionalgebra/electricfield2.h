#ifndef TEMPLAT_LATTICE_ALGEBRA_AXIONALGEBRA_ELECTRICFIELD2_H
#define TEMPLAT_LATTICE_ALGEBRA_AXIONALGEBRA_ELECTRICFIELD2_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"

namespace TempLat
{
  /**
   * @vocab-summary Site-centred electric field: the average of the two links straddling a site in each
   * direction. Returns a vector expression when called without an index.
   * @vocab-signature electricField2(Es)   electricField2(Es, i)
   **/
  template <typename R> auto electricField2(R Es, Tag<1> t) { return 0.5 * (Es(1_c) + shift(Es(1_c), -1_c)); }

  template <typename R> auto electricField2(R Es, Tag<2> t) { return 0.5 * (Es(2_c) + shift(Es(2_c), -2_c)); }

  template <typename R> auto electricField2(R Es, Tag<3> t) { return 0.5 * (Es(3_c) + shift(Es(3_c), -3_c)); }

  template <typename R> auto electricField2(R &&Es) { return MakeVector(i, 1, 3, electricField2(Es, i)); }
} // namespace TempLat

#endif
