#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSSUBTRACT_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSSUBTRACT_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/hasstaticgetter.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessbinaryoperator.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelesswrapper.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/allmatrixcomponents.h"

namespace TempLat
{
  /** @brief A class which implements complex subtraction.
   *
   *
   * Unit test: ctest -R test-complexfieldadd
   **/
  template <class R, class T> class SymTracelessSubtraction : public SymTracelessBinaryOperator<R, T>
  {
  public:
    using SymTracelessBinaryOperator<R, T>::mR;
    using SymTracelessBinaryOperator<R, T>::mT;

    // Put public methods here. These should change very little over time.
    SymTracelessSubtraction(const R &pR, const T &pT) : SymTracelessBinaryOperator<R, T>(pR, pT) {}

    auto SymTracelessGet(Tag<0> t) const { return getComponent(mR, 0_c) - getComponent(mT, 0_c); }
    auto SymTracelessGet(Tag<1> t) const { return getComponent(mR, 1_c) - getComponent(mT, 1_c); }
    auto SymTracelessGet(Tag<2> t) const { return getComponent(mR, 2_c) - getComponent(mT, 2_c); }
    auto SymTracelessGet(Tag<3> t) const { return getComponent(mR, 3_c) - getComponent(mT, 3_c); }
    auto SymTracelessGet(Tag<4> t) const { return getComponent(mR, 4_c) - getComponent(mT, 4_c); }

    auto SymTracelessGet(Tag<1> t1, Tag<1> t2) const { return (*this).SymTracelessGet(0_c); }
    auto SymTracelessGet(Tag<1> t1, Tag<2> t2) const { return (*this).SymTracelessGet(1_c); }
    auto SymTracelessGet(Tag<1> t1, Tag<3> t2) const { return (*this).SymTracelessGet(2_c); }
    auto SymTracelessGet(Tag<2> t1, Tag<1> t2) const { return (*this).SymTracelessGet(1_c); }
    auto SymTracelessGet(Tag<2> t1, Tag<2> t2) const { return (*this).SymTracelessGet(3_c); }
    auto SymTracelessGet(Tag<2> t1, Tag<3> t2) const { return (*this).SymTracelessGet(4_c); }
    auto SymTracelessGet(Tag<3> t1, Tag<1> t2) const { return (*this).SymTracelessGet(2_c); }
    auto SymTracelessGet(Tag<3> t1, Tag<2> t2) const { return (*this).SymTracelessGet(4_c); }
    auto SymTracelessGet(Tag<3> t1, Tag<3> t2) const
    {
      return -(*this).SymTracelessGet(0_c) - (*this).SymTracelessGet(3_c);
    }

    template <int N, int M> auto operator()(Tag<N> t1, Tag<M> t2) const
    {
      static_assert(N >= 1 && N <= 3 && M >= 1 && M <= 3, "Operator(): N and M must be between 1 and 3 for SymWrapper");
      return SymTracelessGet(t1, t2);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, std::decay_t<T> t, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
        DoEval::eval(t, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      auto cL = DoEval::eval(mR, idx...);
      auto cR = DoEval::eval(mT, idx...);
      device::array<decltype(cL[0] - cR[0]), 5> result;
      result[0] = cL[0] - cR[0];
      result[1] = cL[1] - cR[1];
      result[2] = cL[2] - cR[2];
      result[3] = cL[3] - cR[3];
      result[4] = cL[4] - cR[4];
      return result;
    }

    virtual std::string operatorString() const override { return "-"; }
  };

  template <typename R, typename T>
    requires(HasSymTracelessGet<R> && HasSymTracelessGet<T>)
  auto operator-(const R &r, const T &t)
  {
    return SymTracelessSubtraction<R, T>{r, t};
  }

} // namespace TempLat

#endif
