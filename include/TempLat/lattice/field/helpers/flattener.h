#ifndef TEMPLAT_LATTICE_FIELD_HELPERS_FLATTENER_H
#define TEMPLAT_LATTICE_FIELD_HELPERS_FLATTENER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/lattice/field/helpers/hasastuplecat.h"

namespace TempLat
{
  /** @brief A class which flattents a composite collection.
   *
   *
   * Unit test: ctest -R test-flattener
   **/
  class Flattener
  {
  public:
    // Put public methods here. These should change very little over time.

    template <typename U>
      requires HasAsTupleCat<U>::value
    static inline auto asTupleCat(U &obj)
    {
      return obj.asTupleCat();
    }

    template <typename U>
      requires(!HasAsTupleCat<U>::value)
    static inline auto asTupleCat(U &&obj)
    {
      return std::forward<U>(obj);
    }
  };
} // namespace TempLat

#endif
