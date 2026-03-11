#ifndef TEMPLAT_LATTICE_ALGEBRA_COMPLEXALGEBRA_FIELDSTRENGTH_H
#define TEMPLAT_LATTICE_ALGEBRA_COMPLEXALGEBRA_FIELDSTRENGTH_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/util/rangeiteration/tag.h"
#include "TempLat/lattice/algebra/spatialderivatives/forwdiff.h"
#include "TempLat/lattice/algebra/spatialderivatives/neutdiff.h"

namespace TempLat
{
  /** @brief Returns the field strength tensor F_{mu,nu} = d_mu A_nu - d_nu A_mu, where A is a gauge field and d is a
   * forward derivative, using a forward finite difference.
   *
   * @param A The gauge field
   * @param mu The first index of the field strength tensor.
   * @param nu The second index of the field strength tensor.
   **/
  template <typename R, int Mu, int Nu> auto fieldStrength(R A, Tag<Mu> mu, Tag<Nu> nu)
  {
    return forwDiff(A(nu), mu) - forwDiff(A(mu), nu);
  }

  /** @brief Returns the field strength tensor F_{mu,nu} = d_mu A_nu - d_nu A_mu, where A is a gauge field and d is a
   * forward derivative, using a central finite difference.
   *
   * @param A The gauge field
   * @param mu The first index of the field strength tensor.
   * @param nu The second index of the field strength tensor.
   **/
  template <typename R, int Mu, int Nu> auto fieldStrengthCtr(R A, Tag<Mu> mu, Tag<Nu> nu)
  {
    return neutDiff(A(nu), mu) - neutDiff(A(mu), nu);
  }
} // namespace TempLat

#endif
