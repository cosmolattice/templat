#ifndef COSMOINTERFACE_GAUGEALGEBRA_PLAQUETTEBACK_H
#define COSMOINTERFACE_GAUGEALGEBRA_PLAQUETTEBACK_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/su2algebra/su2dagger.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquette.h"
#include "TempLat/lattice/algebra/listoperators/listshift.h"
#include "TempLat/lattice/memory/memorytoolbox.h"

namespace TempLat
{
  /** @brief A class which implemets a backward plaquette. Useful in equations of motions for example.
   *
   *
   * Unit test: ctest -R test-plaquetteback
   *
   * @vocab-summary The plaquette traversed backwards in $\nu$, needed alongside plaq when assembling equations
   * of motion.
   * @vocab-signature plaqBack(Us, mu, nu)
   **/

  template <int Mu, int Nu, typename R> auto plaqBack(const R &Us, Tag<Mu> mu, Tag<Nu> nu)
  {
    return shift(dagger(Us(nu)), -nu) * shift(Us(mu), -nu) * (shift(shift(Us(nu), mu), -nu) * dagger(Us(mu)));
    // return (dagger(shift<-Nu>(Us(nu))) * shift<-Nu>(Us(mu))) * (shift<-Nu>(shift<Mu>(Us(nu)) ) *
    // dagger(shift<Mu>(Us(mu)))) * (dagger(Us(nu))  *shift<-Nu>(Us(nu)));
  }
} // namespace TempLat

#endif
