#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISSCALARTYPE_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_ISSCALARTYPE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler,  Year: 2025

#include "TempLat/parallel/device.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2doubletget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hasmatrixget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymtracelessget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hashermget.h"

namespace TempLat
{
  template <typename T>
  concept IsScalarType =
  (std::is_arithmetic_v<T> || HasEvalMethod<T>) && !HasComplexFieldGet<T> && !HasSU2Get<T> && !HasSU2DoubletGet<T> && !HasSymGet<T> && !HasHermGet<T> && !HasMatrixGet<T> && !HasSymTracelessGet<T>;

} // namespace TempLat

#endif
