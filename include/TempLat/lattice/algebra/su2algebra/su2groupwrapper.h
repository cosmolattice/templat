#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2GROUPWRAPPER_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2GROUPWRAPPER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/lattice/algebra/su2algebra/su2operator.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/lattice/algebra/operators/power.h"

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"

namespace TempLat
{
  /** @brief A class which computes an element of the group SU(2).
   * By that, we mean that det=1 is imposed as a constrained on the 0th element, like in the SU2Fields.
   *
   * Unit test: ctest -R test-su2groupwrapper
   **/
  template <class A, class B, class C> class SU2GroupWrapper : public SU2Operator
  {
  public:
    // Put public methods here. These should change very little over time.
    using SV = typename GetGetReturnType<A>::type;

    SU2GroupWrapper(const A &pA, const B &pB, const C &pC) : mA(pA), mB(pB), mC(pC) {}

    template <int N> auto SU2Get(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3,
                    "SU2Get: N must be greater than 0 and less than or equal to 3 for SU2GroupWrapper");

      if constexpr (N == 0) {
        return sqrt(1.0 - pow<2>(mA) - pow<2>(mB) - pow<2>(mC));
      } else {
        return device::get<N - 1>(device::tie(mA, mB, mC));
      }
    }

    template <typename... IDX>
      requires requires(std::decay_t<A> a, std::decay_t<B> b, std::decay_t<C> c, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(a, idx...);
        DoEval::eval(b, idx...);
        DoEval::eval(c, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<SV, 4> result;
      result[1] = DoEval::eval(mA, idx...);
      result[2] = DoEval::eval(mB, idx...);
      result[3] = DoEval::eval(mC, idx...);
      result[0] = sqrt(SV(1) - powr<2>(result[1]) - powr<2>(result[2]) - powr<2>(result[3]));
      return result;
    }

    template <int N> auto operator()(Tag<N> t) { return SU2Get(t); }

    void preGet()
    {
      PreGet::apply(mA);
      PreGet::apply(mB);
      PreGet::apply(mC);
    }

    void postGet()
    {
      PostGet::apply(mA);
      PostGet::apply(mB);
      PostGet::apply(mC);
    }

    std::string toString() const
    {
      return "SU2Group(" + GetString::get(mA) + "," + GetString::get(mB) + "," + GetString::get(mC) + ")";
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    A mA;
    B mB;
    C mC;
  };

  template <class A, class B, class C> auto SU2GroupWrap(A &&pA, B &&pB, C &&pC)
  {
    return SU2GroupWrapper<A, B, C>(pA, pB, pC);
  }

  template <class R> auto toSU2(R r) { return SU2GroupWrap(r.SU2Get(1_c), r.SU2Get(2_c), r.SU2Get(3_c)); }
} // namespace TempLat

#endif
