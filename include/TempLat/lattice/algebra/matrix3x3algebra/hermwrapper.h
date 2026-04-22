#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HERMWRAPPER_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HERMWRAPPER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/allmatrixoperator.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/complexalgebra/real.h"
#include "TempLat/lattice/algebra/complexalgebra/imag.h"

#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"

namespace TempLat
{
  /** @brief A class which wraps five objects as an hermitian 3x3 matrix  .
   *
   *
   * Unit test: ctest -R test-complexwrapper
   **/
  /*TODO: Jorge: Maybe it would be great to check the diagonal elements are real. For now, I assumme this is checked by
   * the user */
  template <typename R11, typename R12, typename R13, typename R22, typename R23, typename R33>
  class HermWrapper : public HermOperator
  {
  public:
    // Put public methods here. These should change very little over time.

    HermWrapper(const R11 &pR11, const R12 &pR12, const R13 &pR13, const R22 &pR22, const R23 &pR23, const R33 &pR33)
        : mR11(pR11), mR12(pR12), mR13(pR13), mR22(pR22), mR23(pR23), mR33(pR33)
    {
    }

    auto HermGet(Tag<0> t) const { return mR11; }
    auto HermGet(Tag<1> t) const { return mR12; }
    auto HermGet(Tag<2> t) const { return mR13; }
    auto HermGet(Tag<3> t) const { return mR22; }
    auto HermGet(Tag<4> t) const { return mR23; }
    auto HermGet(Tag<5> t) const { return mR33; }

    auto HermGet(Tag<1> t1, Tag<1> t2) const { return mR11; }
    auto HermGet(Tag<1> t1, Tag<2> t2) const { return mR12; }
    auto HermGet(Tag<1> t1, Tag<3> t2) const { return mR13; }
    auto HermGet(Tag<2> t1, Tag<1> t2) const { return conj(mR12); }
    auto HermGet(Tag<2> t1, Tag<2> t2) const { return mR22; }
    auto HermGet(Tag<2> t1, Tag<3> t2) const { return mR23; }
    auto HermGet(Tag<3> t1, Tag<1> t2) const { return conj(mR13); }
    auto HermGet(Tag<3> t1, Tag<2> t2) const { return conj(mR23); }
    auto HermGet(Tag<3> t1, Tag<3> t2) const { return mR33; }

    template <int N> auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 5, "Operator(): N must be between 0 and 5 for HermWrapper");
      return HermGet(t);
    }

    template <int N, int M> auto operator()(Tag<N> t1, Tag<M> t2) const
    {
      static_assert(N >= 1 && N <= 3 && M >= 1 && M <= 3,
                    "Operator(): N and M must be between 1 and 3 for HermWrapper");
      return HermGet(t1, t2);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R11> r11, std::decay_t<R12> r12, std::decay_t<R13> r13, std::decay_t<R22> r22,
                        std::decay_t<R23> r23, std::decay_t<R33> r33, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r11, idx...);
        DoEval::eval(r12, idx...);
        DoEval::eval(r13, idx...);
        DoEval::eval(r22, idx...);
        DoEval::eval(r23, idx...);
        DoEval::eval(r33, idx...);
      }
    DEVICE_INLINE_FUNCTION auto
    eval(const IDX &...idx) const // TODO: Jorge: I need to discuss if the use of aux here takes time or if it is fine.
    {
      device::array<complex<decltype(DoEval::eval(mR11, idx...))>, 6> result;
      result[0] = complex(DoEval::eval(mR11, idx...), static_cast<decltype(DoEval::eval(mR11, idx...))>(0));
      auto aux = DoEval::eval(mR12, idx...);
      result[1] = complex(aux[0], aux[1]);
      aux = DoEval::eval(mR13, idx...);
      result[2] = complex(aux[0], aux[1]);
      result[3] = complex(DoEval::eval(mR22, idx...), static_cast<decltype(DoEval::eval(mR22, idx...))>(0));
      aux = DoEval::eval(mR23, idx...);
      result[4] = complex(aux[0], aux[1]);
      result[5] = complex(DoEval::eval(mR33, idx...), static_cast<decltype(DoEval::eval(mR33, idx...))>(0));
      return result;
    }

    void preGet()
    {
      PreGet::apply(mR11);
      PreGet::apply(mR12);
      PreGet::apply(mR13);
      PreGet::apply(mR22);
      PreGet::apply(mR23);
      PreGet::apply(mR33);
    }

    void postGet()
    {
      PostGet::apply(mR11);
      PostGet::apply(mR12);
      PostGet::apply(mR13);
      PostGet::apply(mR22);
      PostGet::apply(mR23);
      PostGet::apply(mR33);
    }

    std::string toString() const
    {
      return "Herm( (" + GetString::get(mR11) + "," + GetString::get(mR12) + "," + GetString::get(mR13) + ") , (" +
             GetString::get(conj(mR12)) + "," + GetString::get(mR22) + "," + GetString::get(mR23) + ") , (" +
             GetString::get(conj(mR13)) + "," + GetString::get(conj(mR23)) + "," + GetString::get(mR33) + ") ";
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    R11 mR11;
    R12 mR12;
    R13 mR13;
    R22 mR22;
    R23 mR23;
    R33 mR33;
  };

  template <typename R11, typename R12, typename R13, typename R22, typename R23, typename R33>
  HermWrapper<R11, R12, R13, R22, R23, R33> ConstructHerm(const R11 &r11, const R12 &r12, const R13 &r13,
                                                          const R22 &r22, const R23 &r23, const R33 &r33)
  {
    return {r11, r12, r13, r22, r23, r33};
  }
} // namespace TempLat

#endif
