#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_ARG_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_ARG_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/operators/arg.h"
#include "TempLat/lattice/algebra/operators/arg2.h"
#include "TempLat/lattice/algebra/helpers/iscomplextype.h"

namespace TempLat
{
  /** @brief A class which returns teh phase of a complex field, between -pi and pi
   *
   * Unit test: ctest -R test-imag
   **/
  template <class T>
    requires HasComplexFieldGet<T>
  static auto arg2(T &&t)
  {
    return arg2(t.ComplexFieldGet(Tag<1>()), t.ComplexFieldGet(Tag<0>()));
  }

  template <class T>
    requires IsComplexType<T>
  static auto arg2(T &&t)
  {
    return arg2(t.imag(), t.real());
  }

  template <class T>
  requires HasComplexFieldGet<T>
  static auto arg(T &&t)
  {
    return arg(t.ComplexFieldGet(Tag<1>()), t.ComplexFieldGet(Tag<0>()));
  }

  template <class T>
  requires IsComplexType<T>
  static auto arg(T &&t)
  {
    return arg(t.imag(), t.real());
  }
} // namespace TempLat

#endif
