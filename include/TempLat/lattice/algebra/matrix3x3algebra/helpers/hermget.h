#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_HERMGET_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_HERMGET_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/util/rangeiteration/tag.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hashermget.h"

namespace TempLat
{
  /** @brief A class which get real and imaginary part.
   *
   *
   * Unit test: ctest -R test-hermget
   **/
  class HermGetter
  {
    public:
    template <typename R, int N>
      requires(HasHermGet<R>)
    static auto get(R &&r, Tag<N> t)
    {
      return r.HermGet(t);
    }
  };
} // namespace TempLat

#endif
