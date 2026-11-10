#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_GHOSTSHUNTER_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_GHOSTSHUNTER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/helpers/confirmghosts.h"

namespace TempLat
{
  /** @brief Entry point of the ghost walk for a whole expression: checks the ghost depth of every field it reads
   *  at neighbouring sites and updates their ghosts if needed. See ConfirmGhosts.
   *
   * Unit test: ctest -R test-ghostreach
   **/
  class GhostsHunter
  {
  public:
    template <typename U> static void apply(U &obj) { ConfirmGhosts::apply(obj, GhostReach{}); }

  private:
    GhostsHunter() = default;
  };
} // namespace TempLat

#endif
