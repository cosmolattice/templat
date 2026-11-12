#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_CONFIRMGHOSTS_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_CONFIRMGHOSTS_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/helpers/hasghostmethod.h"
#include "TempLat/lattice/algebra/helpers/iscompositeexpression.h"

namespace TempLat
{
  /** @brief The ghost walk. Every node forwards confirmGhostsUpToDate(reach) to its operands; non-local nodes
   *  (shifts, finite differences, Laplacians) first add their stencil to reach. A field reached with a non-zero
   *  reach checks its ghost depth and updates its ghosts; see GhostReach. Objects without the method (constants)
   *  are skipped. Returns the number of FFTs and ghost exchanges triggered.
   *
   * Unit test: ctest -R test-ghostreach
   **/
  class ConfirmGhosts
  {
  public:
    template <typename U>
      requires HasGhostMethod<U>
    static device::Idx apply(U &obj, const GhostReach &reach)
    {
      return obj.confirmGhostsUpToDate(reach);
    }

    template <typename U>
      requires(!HasGhostMethod<U>)
    static constexpr device::Idx apply(U &, const GhostReach &)
    {
      static_assert(!IsCompositeExpression<U>,
                    "Composite expression node without confirmGhostsUpToDate: its operands would skip the ghost "
                    "exchange.");
      return 0;
    }

    template <int N, typename U>
      requires HasGhostMethodIndexed<N, U>
    static device::Idx apply(U &obj, Tag<N> i, const GhostReach &reach)
    {
      return obj.confirmGhostsUpToDate(i, reach);
    }

    template <int N, typename U>
      requires(!HasGhostMethodIndexed<N, U> && HasGhostMethodElement<N, U>)
    static device::Idx apply(U &obj, Tag<N> i, const GhostReach &reach)
    {
      return obj(i).confirmGhostsUpToDate(reach);
    }

    template <int N, typename U>
      requires(!HasGhostMethodIndexed<N, U> && !HasGhostMethodElement<N, U>)
    static constexpr device::Idx apply(U &, Tag<N>, const GhostReach &)
    {
      return 0;
    }

  private:
    ConfirmGhosts() = default;
  };
} // namespace TempLat

#endif
