#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_MATRIXGETGETRETURNTYPE_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_MATRIXGETGETRETURNTYPE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026

#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{
  /** @brief A class which unpack the return type of a complex type.
   *
   * Unit test: ctest -R test-complexgetgetreturntype
   **/
  template <typename T> struct SymTracelessGetGetReturnType {
    using type = typename GetGetReturnType<
        std::decay_t<decltype(std::declval<T>().SymTracelessGet(std::declval<Tag<0>>()))>>::type;
  };

  template <typename T> struct SymGetGetReturnType {
    using type =
        typename GetGetReturnType<std::decay_t<decltype(std::declval<T>().SymGet(std::declval<Tag<0>>()))>>::type;
  };

  template <typename T> struct HermGetGetReturnType {
    using type =
        typename GetGetReturnType<std::decay_t<decltype(std::declval<T>().HermFieldGet(std::declval<Tag<0>>()))>>::type;
  };

  template <typename T> struct MatrixGetGetReturnType {
    using type = typename GetGetReturnType<
        std::decay_t<decltype(std::declval<T>().MatrixFieldGet(std::declval<Tag<0>>()))>>::type;
  };

} // namespace TempLat

#endif
