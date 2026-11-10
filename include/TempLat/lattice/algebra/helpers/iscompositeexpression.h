#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISCOMPOSITEEXPRESSION_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISCOMPOSITEEXPRESSION_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <type_traits>

#include "TempLat/lattice/algebra/helpers/ghostreach.h"
#include "TempLat/parallel/device.h"

#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2doubletget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymtracelessget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hashermget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hasmatrixget.h"

namespace TempLat
{
  /** @brief An expression of one of the composite algebras (complex, SU(2), SU(2) doublet, 3x3 matrices).
   *
   * Assignments walk these expressions as a whole for space and ghost confirmation, so every composite node must
   * forward confirmGhostsUpToDate and confirmSpace to its operands. A node that does not would silently skip the
   * ghost exchange of everything below it; ConfirmGhosts and ConfirmSpace therefore reject composite nodes without
   * these methods at compile time. Constants derive from NoFieldDependence.
   *
   * Unit test: ctest -R test-assignmentaliasing
   **/
  template <class T>
  concept IsCompositeExpression =
      HasComplexFieldGet<std::decay_t<T>> || HasSU2Get<std::decay_t<T>> || HasSU2DoubletGet<std::decay_t<T>> ||
      HasSymTracelessGet<std::decay_t<T>> || HasSymGet<std::decay_t<T>> || HasHermGet<std::decay_t<T>> ||
      HasMatrixGet<std::decay_t<T>>;

  /** @brief Base for constant composite expressions: they read no field, so there is nothing to confirm. */
  struct NoFieldDependence {
    static device::Idx confirmGhostsUpToDate(const GhostReach &) { return 0; }
    template <typename Layout, typename SpaceType> static void confirmSpace(const Layout &, const SpaceType &) {}
  };
} // namespace TempLat

#endif
