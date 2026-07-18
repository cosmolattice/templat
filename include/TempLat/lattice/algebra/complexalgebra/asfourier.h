#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_ASFOURIER_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_ASFOURIER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/lattice/algebra/operators/unaryoperator.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/complexgetgetreturntype.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

namespace TempLat
{
  /** @brief A class which treats a complex field as an object in fourier space.
   *
   * Derives from UnaryOperator so that the pre/postGet, confirmSpace, ghost and toolBox lifecycle is
   * forwarded to the wrapped child. Without this the wrapped subtree is silently skipped (e.g. a random
   * field never advances its generation), so eval alone is not enough.
   *
   * Unit test: ctest -R test-asfourier
   **/
  template <typename R> class ComplexFieldAsFourier : public UnaryOperator<R>
  {
  public:
    using UnaryOperator<R>::mR;

    // Put public methods here. These should change very little over time.
    using mRType = typename ComplexGetGetReturnType<R>::type;

    ComplexFieldAsFourier(const R &pR) : UnaryOperator<R>(pR) {}

    static constexpr size_t NDim = GetNDim::get<R>();

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto result = DoEval::eval(mR, idx...);
      return complex<mRType>(result[0], result[1]);
    }

    std::string toString() const { return "asFourier(" + GetString::get(mR) + ")"; }
  };

  template <typename R> ComplexFieldAsFourier<R> asFourier(R &&r)
  {
    return ComplexFieldAsFourier<R>(std::forward<R>(r));
  }
} // namespace TempLat

#endif
