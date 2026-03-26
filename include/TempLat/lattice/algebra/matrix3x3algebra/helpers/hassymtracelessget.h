#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_HASSYMTRACELESSFIELDGET_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HELPERS_HASSYMTRACELESSFIELDGET_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{

  /** @brief A concept which checks if instance is a symmetric traceless matrix.
   *
   * Unit test: ctest -R test-hassymmetrictracelessfieldget
   **/
  template <class T>
  concept HasSymTracelessGet1 = requires(T t, Tag<0> tag) { t.SymTracelessGet(tag); };

  template <class T>
  concept HasSymTracelessGet2 = requires(T t, Tag<1> tag1, Tag<1> tag2) { t.SymTracelessGet(tag1, tag2); };

  template <class T>
  concept HasSymTracelessGet = HasSymTracelessGet1<T> || HasSymTracelessGet2<T>;
} // namespace TempLat

#endif
