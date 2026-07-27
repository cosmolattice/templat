#ifndef COSMOINTERFACE_SU2ALGEBRA_SU2DAGGER_H
#define COSMOINTERFACE_SU2ALGEBRA_SU2DAGGER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2019

#include <type_traits>

#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/su2algebra/su2unaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2multiply.h"
#include "TempLat/lattice/algebra/su2algebra/su2shift.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"

namespace TempLat
{
  /** @brief A class which implements the hermitian conjugation.
   *
   *
   * Unit test: ctest -R test-su2dagger
   **/
  template <typename R> class SU2Dagger : public SU2UnaryOperator<R>
  {
  public:
    using SV = typename SU2GetGetReturnType<R>::type;
    using SU2UnaryOperator<R>::mR;

    // Put public methods here. These should change very little over time.
    SU2Dagger(const R &pR) : SU2UnaryOperator<R>(pR) {}

    template <int N> auto SU2Get(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "SU2Get: N must be between 0 and 3 for SU2Dagger");
      if constexpr (N == 0) {
        return mR.SU2Get(Tag<N>());
      } else {
        return -mR.SU2Get(Tag<N>());
      }
    }

    template <int N> auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "Operator(): N must be between 0 and 3 for SU2Dagger");
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
      c[1] = -c[1];
      c[2] = -c[2];
      c[3] = -c[3];
      return c;
    }

    std::string toString() const { return GetString::get(mR) + "^\u2020"; }
  };

  /**
   * @vocab-summary Hermitian conjugate. For SU(2) it negates the three vector components; the operator factory
   * also simplifies as it builds, collapsing $(A^\dagger)^\dagger$ to $A$ and commuting past a shift.
   **/
  template <class R>
    requires HasSU2Get<R>
  auto dagger(const R &r)
  {
    return SU2Dagger<R>(r);
  };

  // ---------------------------------------------------------------------------------------------------
  // Algebraic simplification of dagger expressions.
  //
  // SU(2) dagger is the quaternion conjugate, an anti-involution: (U^dagger)^dagger = U and
  // (A.B)^dagger = B^dagger.A^dagger, valid for any quaternion-valued expression (unit or scaled). We
  // apply these ONLY when they remove a dagger (i.e. when a sub-dagger cancels): the plain
  // anti-homomorphism dag(A.B) -> dag(B).dag(A) on two non-dagger operands would turn one 3-component
  // negation into two, a pessimisation, so we never fire it there. The dag(dag U) -> U case is also
  // folded by the optimiser at runtime, but eliminating it here additionally shrinks the instantiated
  // type. Net effect: removes daggers the back-end optimiser cannot (it does not know the quaternion
  // anti-homomorphism), and exposes leaf-level operands for common-subexpression elimination.
  // ---------------------------------------------------------------------------------------------------

  template <class T> struct IsSU2DaggerT : std::false_type {
  };
  template <class R> struct IsSU2DaggerT<SU2Dagger<R>> : std::true_type {
  };
  template <class T>
  concept IsSU2Dagger = IsSU2DaggerT<std::decay_t<T>>::value;

  // Involution: (U^dagger)^dagger = U.
  template <class R> auto dagger(const SU2Dagger<R> &r) { return r.getOperand(); }

  // Anti-homomorphism, fired only when at least one operand is itself a dagger, so that the rewrite
  // cancels that dagger through the involution above: (A.B)^dagger = B^dagger.A^dagger.
  template <class A, class B>
    requires(IsSU2Dagger<A> || IsSU2Dagger<B>)
  auto dagger(const SU2Multiplication<A, B> &m)
  {
    return dagger(m.getSecond()) * dagger(m.getFirst());
  }

  // Commute dagger past a shift: (shift(X))^dagger = shift(X^dagger). A shift is a coordinate relabeling
  // and the dagger a per-site quaternion conjugation, so they commute; cost-neutral at eval. Unlike the
  // anti-homomorphism this fires unconditionally - it lets nested dag(shift(dag(Y))) collapse to
  // shift(Y) via the involution above (e.g. inside the axion-current theta_ijk).
  template <class R, int... N> auto dagger(const SU2Shifter<R, N...> &s) { return shift<N...>(dagger(s.getOperand())); }
  template <class R, int N> auto dagger(const SU2ShifterByOne<R, N> &s) { return shift<N>(dagger(s.getOperand())); }

  /**
   * @vocab-summary Short spelling of dagger.
   **/
  template <class R>
    requires HasSU2Get<R>
  auto dag(const R &r)
  {
    return dagger(r);
  };
} // namespace TempLat

#endif
