#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_HASGHOSTMETHOD_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_HASGHOSTMETHOD_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2026

#include "TempLat/util/rangeiteration/tag.h"

#include "TempLat/lattice/algebra/helpers/ghostreach.h"
namespace TempLat
{
  /** @brief Concepts for the ghost walk: whether an object takes part in it via confirmGhostsUpToDate(reach), as a
   *  whole, per list component (Tag<N>), or through its N-th element.
   *
   * Unit test: ctest -R test-hasghostmethod
   **/
  template <class T>
  concept HasGhostMethod = requires(std::decay_t<T> t, const GhostReach &reach) { t.confirmGhostsUpToDate(reach); };

  template <int N, class T>
  concept HasGhostMethodIndexed =
      requires(std::decay_t<T> t, Tag<N> tag, const GhostReach &reach) { t.confirmGhostsUpToDate(tag, reach); };

  template <int N, class T>
  concept HasGhostMethodElement =
      requires(std::decay_t<T> t, Tag<N> tag, const GhostReach &reach) { t(tag).confirmGhostsUpToDate(reach); };
} // namespace TempLat

#endif
