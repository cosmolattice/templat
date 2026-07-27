#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_IMAG_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_IMAG_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/complexalgebra/ascomplexfield.h"
#include "TempLat/lattice/algebra/helpers/iscomplextype.h"

namespace TempLat
{
  /** @brief A class which returns imaginary part of a fields or complex number.
   *
   * Unit test: ctest -R test-imag
   *
   * @vocab-summary Imaginary part of a complex field or complex number.
   **/
  template <class T>
    requires HasComplexFieldGet<T>
  auto Imag(T &&t)
  {
    return t.ComplexFieldGet(Tag<1>());
  }

  template <class T>
    requires IsComplexType<T>
  auto Imag(T &&t)
  {
    return t.imag();
  }

  /** @brief Imaginary part of a getter whose eval() returns a complex value. Bridges it into the
   *  complex-field protocol first, so Imag(f) == Imag(asComplexField(f)). */
  template <class T>
    requires HasComplexEval<std::decay_t<T>>
  auto Imag(T &&t)
  {
    return asComplexField(std::forward<T>(t)).ComplexFieldGet(1_c);
  }
} // namespace TempLat

#endif
