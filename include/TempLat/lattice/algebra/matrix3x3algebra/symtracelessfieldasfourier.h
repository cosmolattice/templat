#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSFIELDASFOURIER_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSFIELDASFOURIER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026

#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/matrixgetgetreturntype.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

namespace TempLat
{
  /** @brief A class which treats a symmetric-traceless field as an object in fourier space.
   *
   * Unit test: ctest -R test-asfourier
   **/
  template <typename R> class SymTracelessFieldAsFourier
  {
  public:
    // Put public methods here. These should change very little over time.
    using mRType = typename SymTracelessGetGetReturnType<R>::type;

    SymTracelessFieldAsFourier(const R &pR) : mR(pR) {}

    static constexpr size_t NDim = GetNDim::get<R>();

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto result = DoEval::eval(mR, idx...);
      return complex<mRType>(result[0], result[1], result[2], result[3], result[4]);
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    R mR;
  };

  /**
   * @vocab-summary Views a symmetric-traceless field in Fourier space.
   * @vocab-tags SymTraceless, Fourier
   **/
  template <typename R> SymTracelessFieldAsFourier<R> symTracelessFieldAsFourier(R &&r)
  {
    return SymTracelessFieldAsFourier<R>(std::forward<R>(r));
  }
} // namespace TempLat

#endif
