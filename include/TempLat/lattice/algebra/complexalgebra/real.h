#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_REAL_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_REAL_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/complexalgebra/ascomplexfield.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/helpers/iscomplextype.h"

namespace TempLat
{
  /** @brief A class which get real parts of fields.
   *
   *
   * Unit test: ctest -R test-real
   *
   * @vocab-summary Real part of a complex field or complex number.
   **/
  template <class T>
    requires HasComplexFieldGet<T>
  auto Real(T &&t)
  {
    return t.ComplexFieldGet(0_c);
  }

  template <class T>
    requires IsComplexType<T>
  auto Real(T &&t)
  {
    return t.real();
  }

  /** @brief Real part of a getter whose eval() returns a complex value (e.g. RandomGaussianField in
   *  Fourier space). Bridges it into the complex-field protocol first, so Real(f) == Real(asComplexField(f)). */
  template <class T>
    requires HasComplexEval<std::decay_t<T>>
  auto Real(T &&t)
  {
    return asComplexField(std::forward<T>(t)).ComplexFieldGet(0_c);
  }
} // namespace TempLat

#endif
