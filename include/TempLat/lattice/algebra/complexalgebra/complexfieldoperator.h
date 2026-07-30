#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELDOPERATOR_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELDOPERATOR_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/algebra/complexalgebra/helpers/complexfieldget.h"

namespace TempLat
{
  /** @brief A class which implements common features of complex fields operators.
   *
   *
   * Unit test: ctest -R test-complexfieldoperator
   **/
  class ComplexFieldOperator
  {
  public:
    // Put public methods here. These should change very little over time.
    ComplexFieldOperator() = default;

    static constexpr size_t size = 2;
    using Getter = ComplexFieldGetter;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };
} // namespace TempLat

#endif
