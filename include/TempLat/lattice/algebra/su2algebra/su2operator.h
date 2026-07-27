#ifndef COSMOINTERFACE_SU2ALGEBRA_SU2OPERATOR_H
#define COSMOINTERFACE_SU2ALGEBRA_SU2OPERATOR_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/su2algebra/helpers/su2get.h"

namespace TempLat
{
  /** @brief A class which contains basic info about su2 algebra.
   *
   * Unit test: ctest -R test-su2operator
   **/
  class SU2Operator
  {
  public:
    // Put public methods here. These should change very little over time.
    SU2Operator() = default;

    static constexpr size_t size = 4;
    using Getter = SU2Getter;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };
} // namespace TempLat

#endif
