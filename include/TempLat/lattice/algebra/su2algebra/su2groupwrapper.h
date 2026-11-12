#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2GROUPWRAPPER_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2GROUPWRAPPER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/lattice/algebra/su2algebra/su2operator.h"
#include "TempLat/lattice/algebra/su2algebra/su2unaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/paulivectorsalgebra.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/lattice/algebra/operators/power.h"
#include "TempLat/lattice/algebra/operators/squareroot.h"
#include "TempLat/lattice/algebra/operators/divide.h"

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/confirmghosts.h"

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
    SU2GroupWrapper() = default;

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

    /** @brief Space/ghost confirmation forwarded to the wrapped expressions. */
    device::Idx confirmGhostsUpToDate(const GhostReach &reach) const
    {
      return ConfirmGhosts::apply(mA, reach) + ConfirmGhosts::apply(mB, reach) + ConfirmGhosts::apply(mC, reach);
    }
    template <size_t NDim> void confirmSpace(const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      ConfirmSpace::apply(mA, newLayout, spaceType);
      ConfirmSpace::apply(mB, newLayout, spaceType);
      ConfirmSpace::apply(mC, newLayout, spaceType);
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

  /**
   * @vocab-summary Builds a genuine group element from three components, recovering $c_0$ from the unitarity
   * constraint $c_0^2+c_1^2+c_2^2+c_3^2=1$ so that $\det U = 1$ holds by construction. Always gives $c_0 \geq 0$;
   * to keep the sign of an existing element use toSU2.
   * @vocab-signature SU2GroupWrap(c1, c2, c3)
   **/
  template <class A, class B, class C> auto SU2GroupWrap(A &&pA, B &&pB, C &&pC)
  {
    return SU2GroupWrapper<A, B, C>(pA, pB, pC);
  }

  /** @brief Projects an SU(2)-valued expression onto the group, c_mu -> c_mu / |c|.
   *
   * Unit test: ctest -R test-su2groupwrapper
   **/
  template <typename R> class SU2Projector : public SU2UnaryOperator<R>
  {
  public:
    using SU2UnaryOperator<R>::mR;

    SU2Projector(const R &pR) : SU2UnaryOperator<R>(pR) {}

    template <int N> auto SU2Get(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "SU2Get: N must be between 0 and 3 for SU2Projector");
      return mR.SU2Get(t) / sqrt(pow<2>(mR.SU2Get(0_c)) + pow<2>(mR.SU2Get(1_c)) + pow<2>(mR.SU2Get(2_c)) +
                                 pow<2>(mR.SU2Get(3_c)));
    }

    template <int N> auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "Operator(): N must be between 0 and 3 for SU2Projector");
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
      PauliVectorsAlgebra::normalize_inplace(c);
      return c;
    }

    std::string toString() const { return "toSU2(" + GetString::get(mR) + ")"; }
  };

  /**
   * @vocab-summary Projects an SU(2)-valued expression onto the group by normalizing all four components,
   *   $c_\mu \to c_\mu / |c|$. Keeps the sign of $c_0$.
   **/
  template <class R>
    requires HasSU2Get<R>
  auto toSU2(const R &r)
  {
    return SU2Projector<R>(r);
  }
} // namespace TempLat

#endif
