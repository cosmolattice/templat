#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_COMPLEXGETGETRETURNTYPE_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_HELPERS_COMPLEXGETGETRETURNTYPE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{
  /** @brief A class which unpack the return type of a complex type.
   *
   * Unit test: ctest -R test-complexgetgetreturntype
   **/
  template <typename T> struct ComplexGetGetReturnType {
    using type = typename GetGetReturnType<
        std::decay_t<decltype(std::declval<T>().ComplexFieldGet(std::declval<Tag<0>>()))>>::type;

    static constexpr bool isComplex = IsComplexType<type>;
  };

} // namespace TempLat

#endif
