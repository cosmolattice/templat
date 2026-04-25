#ifndef TEMPLAT_LATTICE_LATTICEBCSPEC_H
#define TEMPLAT_LATTICE_LATTICEBCSPEC_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/lattice/ghostcells/boundaryconditions.h"

namespace TempLat
{

  /** @brief Sibling-of-LatticeParameters carrier for a per-dimension boundary-condition spec.
   *
   * Kept separate from `LatticeParameters<T>` so the T-only parameter type stays free of an extra
   * template parameter (which would ripple into every `Field` / `ConfigView` / `AbstractField`
   * ctor signature). User code parses a BC string from config into a `LatticeBCSpec<NDim>` and
   * passes `lBC.bc` to the `Field` / `SU2Field` ctor.
   *
   * Unit test: ctest -R test-latticebcspec
   **/
  template <size_t NDim> struct LatticeBCSpec {
    BCSpec<NDim> bc;

    LatticeBCSpec() : bc(allPeriodic<NDim>()) {}
    explicit LatticeBCSpec(BCSpec<NDim> pBc) : bc(pBc) {}
  };
} // namespace TempLat

#endif
