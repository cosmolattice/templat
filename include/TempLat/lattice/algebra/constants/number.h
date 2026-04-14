#ifndef TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_NUMBER_H
#define TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_NUMBER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/measuringtools/averager.h"

#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A runtime-mutable scalar value that participates in the expression template algebra.
   *
   * Number<T> reports GetNDim = 0 (no NDim member, no getNDim()), causing all spatial
   * derivatives to return ZeroType at the free-function level.
   *
   * It satisfies HasEvalMethod (works in binary expressions) and IsScalarType
   * (no complex/SU2/matrix getters), so Number<T> * Field dispatches correctly.
   *
   * The 4-way operator+= dispatch handles:
   *   1. Lattice expressions (NDim > 0): spatially average, then add
   *   2. 0-dim expressions (another Number, scalar*Number): eval directly
   *   3. Arithmetic (plain double/float): add directly
   *   4. ZeroType: no-op
   *
   * Unit test: ctest -R test-number
   **/
  template <typename T>
  struct Number {
    T value;

    template <typename... IDX>
      requires IsVariadicIndex<IDX...>
    DEVICE_FORCEINLINE_FUNCTION constexpr auto eval(const IDX &...) const
    {
      return value;
    }

    // 1. Lattice expression (NDim > 0): spatially average, then add
    template <class Expr>
      requires(HasEvalMethod<std::decay_t<Expr>> && GetNDim::get<std::decay_t<Expr>>() > 0)
    Number &operator+=(Expr &&e)
    {
      value += average(std::forward<Expr>(e));
      return *this;
    }

    // 2. 0-dim expression (another Number, scalar*Number, etc): eval directly
    template <class Expr>
      requires(HasEvalMethod<std::decay_t<Expr>> && GetNDim::get<std::decay_t<Expr>>() == 0 &&
               !std::is_arithmetic_v<std::decay_t<Expr>>)
    Number &operator+=(Expr &&e)
    {
      value += DoEval::eval(std::forward<Expr>(e), size_t{0});
      return *this;
    }

    // 3. Arithmetic (plain double/float)
    Number &operator+=(T v) { value += v; return *this; }

    // 4. ZeroType: no-op
    Number &operator+=(ZeroType) { return *this; }

    Number &operator=(T v) { value = v; return *this; }
  };
} // namespace TempLat

#endif
