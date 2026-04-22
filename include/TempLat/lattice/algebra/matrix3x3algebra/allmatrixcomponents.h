#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_MATRIXCOMPONENTS_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_MATRIXCOMPONENTS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymtracelessget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hashermget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hasmatrixget.h"
#include "TempLat/lattice/algebra/helpers/iscomplextype.h"
#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{
  /** @brief A class which returns imaginary part of a fields or symmetric-traceless number.
   *
   * Unit test: ctest -R test-imag
   **/
  template <class T, int I>
    requires HasMatrixGet<T>
  static auto getComponent(T &&t, Tag<I> i)
  {
    return t.MatrixGet(i);
  }

  template <class T, int I>
    requires HasSymTracelessGet<T>
  static auto getComponent(T &&t, Tag<I> i)
  {
    return t.SymTracelessGet(i);
  }

  template <class T, int I>
    requires HasSymGet<T>
  static auto getComponent(T &&t, Tag<I> i)
  {
    return t.SymGet(i);
  }

  template <class T, int I>
    requires HasHermGet<T>
  static auto getComponent(T &&t, Tag<I> i)
  {
    return t.HermGet(i);
  }

  template <class T, int I, int J>
    requires HasMatrixGet<T>
  static auto getComponent(T &&t, Tag<I> i, Tag<J> j)
  {
    return t.MatrixGet(i, j);
  }

  template <class T, int I, int J>
    requires HasHermGet<T>
  static auto getComponent(T &&t, Tag<I> i, Tag<J> j)
  {
    return t.HermGet(i, j);
  }

  template <class T, int I, int J>
    requires HasSymGet<T>
  static auto getComponent(T &&t, Tag<I> i, Tag<J> j)
  {
    return t.SymGet(i, j);
  }

  template <class T, int I, int J>
    requires HasSymTracelessGet<T>
  static auto getComponent(T &&t, Tag<I> i, Tag<J> j)
  {
    return t.SymTracelessGet(i, j);
  }

} // namespace TempLat

#endif
