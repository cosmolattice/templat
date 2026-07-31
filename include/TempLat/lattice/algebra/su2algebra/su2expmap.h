#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2EXPMAP_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2EXPMAP_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2024

#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/paulivectorsalgebra.h"
#include "TempLat/lattice/algebra/operators/cosine.h"
#include "TempLat/lattice/algebra/operators/sine.h"
#include "TempLat/lattice/algebra/operators/squareroot.h"
#include "TempLat/lattice/algebra/operators/power.h"
#include "TempLat/lattice/algebra/operators/divide.h"
#include "TempLat/lattice/algebra/su2algebra/su2unaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"

#include "TempLat/parallel/device.h"

namespace TempLat
{

  /** @brief A class which computes the exponential map for su(2).
   *
   *
   * Unit test: ctest -R test-su2expmap
   **/

  template <typename R> class SU2ExpMap : public SU2UnaryOperator<R>
  {
  public:
    using SU2UnaryOperator<R>::mR;

    SU2ExpMap(const R &pR) : SU2UnaryOperator<R>(pR) {}

    template <int N> auto SU2Get(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "SU2Get: N must be between 0 and 3 for SU2ExpMap");

      if constexpr (N == 0) {
        return cos(sqrt(pow<2>(mR.SU2Get(1_c)) + pow<2>(mR.SU2Get(2_c)) + pow<2>(mR.SU2Get(3_c))));
      } else {
        const auto a = sqrt(pow<2>(mR.SU2Get(1_c)) + pow<2>(mR.SU2Get(2_c)) + pow<2>(mR.SU2Get(3_c)));
        return mR.SU2Get(t) / a * sin(a);
      }
    }

    template <int N> auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "Operator(): N must be between 0 and 3 for SU2ExpMap");
      return SU2Get(t);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      // Keep the element type of the evaluated operand; see the note in su2multiply.h
      auto cL = DoEval::eval(mR, idx...);
      PauliVectorsAlgebra::expmap_inplace(cL);
      return cL;
    }

    std::string toString() const { return "exp(" + GetString::get(mR) + ")"; }
  };

  /**
   * @vocab-summary Exponential map from the su(2) algebra to the SU(2) group.
   **/
  template <class R>
    requires HasSU2Get<R>
  auto exp(const R &r)
  {
    return SU2ExpMap<R>(r);
  };

} // namespace TempLat

#endif
