#ifndef COSMOINTERFACE_SU2ALGEBRA_HELPERS_SU2DOUBLETGET_H
#define COSMOINTERFACE_SU2ALGEBRA_HELPERS_SU2DOUBLETGET_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{
  /** @brief A class which return the SU2DouletGet method.
   *
   * Unit test: ctest -R test-su2doubletget
   **/
  class SU2DoubletGetter
  {
  public:
    template <typename R, int N> static auto get(R &&r, Tag<N> t) { return r.SU2DoubletGet(t); }
  };

} // namespace TempLat

#endif
