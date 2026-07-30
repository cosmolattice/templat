#ifndef COSMOINTERFACE_SU2ALGEBRA_HELPERS_SU2DOUBLETGETGETRETURNTYPE_H
#define COSMOINTERFACE_SU2ALGEBRA_HELPERS_SU2DOUBLETGETGETRETURNTYPE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/util/rangeiteration/tag.h"
#include "su2doubletget.h"

namespace TempLat
{
  /** @brief A class which gives the return type of the get function of one of the doublet element.
   *
   * Unit test: ctest -R test-su2doubletgetgetreturntype
   **/
  template <typename T> struct SU2DoubletGetGetReturnType {
    using type = typename GetGetReturnType<
        std::decay_t<decltype(SU2DoubletGetter::get(std::declval<T>(), std::declval<Tag<0>>()))>>::type;
    static constexpr bool isComplex = IsComplexType<type>;
  };
} // namespace TempLat

#endif
