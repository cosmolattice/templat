#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_GETDX_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_GETDX_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/helpers/hasdx.h"
#include "TempLat/parallel/device.h"
namespace TempLat
{
  /** @brief A getter for dx.
   *
   * Unit test: ctest -R test-getdx
   **/
  class GetDx
  {
  public:
    // Put public methods here. These should change very little over time.

    template <typename U>
      requires HasDx<U>
    static auto getDx(U &&obj)
    {
      return obj.getDx();
    }

    template <typename U>
      requires(!HasDx<U>)
    static constexpr int getDx(U &obj)
    {
      return 1;
    }
  };
} // namespace TempLat

#endif
