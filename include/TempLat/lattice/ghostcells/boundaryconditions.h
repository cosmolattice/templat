#ifndef TEMPLAT_LATTICE_GHOSTCELLS_BOUNDARYCONDITIONS_H
#define TEMPLAT_LATTICE_GHOSTCELLS_BOUNDARYCONDITIONS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include <array>
#include <cstddef>

namespace TempLat
{
  enum class BCType { Periodic, Antiperiodic, Dirichlet, Neumann };

  template <size_t NDim> using BCSpec = std::array<BCType, NDim>;

  template <size_t NDim> constexpr BCSpec<NDim> allPeriodic()
  {
    BCSpec<NDim> s{};
    for (auto &b : s) b = BCType::Periodic;
    return s;
  }
} // namespace TempLat

#endif
