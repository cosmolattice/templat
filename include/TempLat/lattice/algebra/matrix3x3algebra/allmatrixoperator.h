#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_MATRIXOPERATOR_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_MATRIXOPERATOR_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/matrixget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/symtracelessget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/symget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hermget.h"

namespace TempLat
{
  /** @brief A class which implements common features of matrix  fields operators.
   *
   *
   * Unit test: ctest -R test-matrixfieldoperator
   **/
  class MatrixOperator
  {
  public:
    // Put public methods here. These should change very little over time.
    DEVICE_FUNCTION
    MatrixOperator() = default;
    static constexpr size_t size = 9;
    using Getter = MatrixGetter;
  };

  class SymOperator
  {
  public:
    // Put public methods here. These should change very little over time.
    DEVICE_FUNCTION
    SymOperator() = default;
    static constexpr size_t size = 6;
    using Getter = SymGetter;
  };

  class HermOperator
  {
  public:
    // Put public methods here. These should change very little over time.
    DEVICE_FUNCTION
    HermOperator() = default;
    static constexpr size_t size = 6;
    using Getter = HermGetter;
  };

  class SymTracelessOperator
  {
  public:
    // Put public methods here. These should change very little over time.
    DEVICE_FUNCTION
    SymTracelessOperator() = default;
    static constexpr size_t size = 5;
    using Getter = SymTracelessGetter;
  };
} // namespace TempLat

#endif
