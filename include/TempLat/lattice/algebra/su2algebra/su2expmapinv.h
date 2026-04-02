#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2EXPMAPINV_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2EXPMAPINV_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2024

#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/operators/power.h"
#include "TempLat/lattice/algebra/operators/cosine.h"
#include "TempLat/lattice/algebra/operators/divide.h"
#include "TempLat/lattice/algebra/operators/acos.h"
#include "TempLat/lattice/algebra/operators/sine.h"
#include "TempLat/lattice/algebra/operators/squareroot.h"
#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/algebra/su2algebra/su2unaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"

#include "TempLat/parallel/device.h"

namespace TempLat
{

  /** @brief A class which computes the inverse exponential map for SU2 (log).
   *
   *
   * Unit test: ctest -R test-su2expmapinv
   **/
  template <typename R> class SU2ExpMapInv : public SU2UnaryOperator<R>
  {
  public:
    using SV = typename SU2GetGetReturnType<R>::type;
    using SU2UnaryOperator<R>::mR;

    SU2ExpMapInv(const R &pR) : SU2UnaryOperator<R>(pR) {}

    template <int N> DEVICE_INLINE_FUNCTION auto SU2Get(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "SU2Get: N must be between 0 and 3 for SU2ExpMapInv");

      if constexpr (N == 0) {
        return ZeroType();
      } else {
        const auto a = acos(mR.SU2Get(0_c));
        return mR.SU2Get(t) * a / sin(a);
      }
    }

    template <int N> DEVICE_INLINE_FUNCTION auto SU2LieAlgebraGet(Tag<N> t) const
    {
      static_assert(N >= 1 && N <= 3, "SU2LieAlgebraGet: N must be between 1 and 3 for SU2ExpMapInv");
      return 2 * SU2Get(t);
    }

    template <int N> DEVICE_INLINE_FUNCTION auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "Operator(): N must be between 0 and 3 for SU2ExpMapInv");
      return SU2Get(t);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      auto c = DoEval::eval(mR, idx...);

      const auto a = device::acos(c[0]);
      const auto sina = device::sin(a);
      // Safe-divide guard: a/sin(a) → 1 as a → 0
      const auto ratio = (sina * sina > SV(1e-30)) ? a / sina : SV(0);

      // We can work in-place on c.
      c[0] = SV(0);
      c[1] = c[1] * ratio;
      c[2] = c[2] * ratio;
      c[3] = c[3] * ratio;
      return c;
    }

    std::string toString() const { return "expinv(" + GetString::get(mR) + ")"; }

    static constexpr size_t numberToSkipAsTuple = 1;
  };

  template <class R>
    requires HasSU2Get<R>
  auto expinv(const R &r)
  {
    return SU2ExpMapInv<R>(r);
  };

} // namespace TempLat

#endif
