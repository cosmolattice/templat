#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_MATRIXGET_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_MATRIXGET_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026

#include "TempLat/util/rangeiteration/tag.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hasmatrixget.h"

namespace TempLat
{
  /** @brief A class which get real and imaginary part.
   *
   *
   * Unit test: ctest -R test-matrixget
   **/
  class MatrixGetter
  {
  public:
    template <typename R, int N>
      requires(HasMatrixGet<R>)
    static auto get(R &&r, Tag<N> t)
    {
      return r.MatrixGet(t);
    }
  };
} // namespace TempLat

#endif
