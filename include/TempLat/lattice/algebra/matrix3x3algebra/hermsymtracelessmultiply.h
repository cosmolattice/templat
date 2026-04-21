#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_HERMSYMTRACELESSMULTIPLY_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_HERMSYMTRACELESSMULTIPLY_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessbinaryoperator.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hermget.h"
#include "TempLat/lattice/algebra/operators/power.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/subtract.h"

namespace TempLat
{
  /** @brief A class which multiplies two symmetric-traceless fields.
   *
   *
   * Unit test: ctest -R test-symmetrictracelessfieldmultiply
   **/
  template <class R, class T> class HermSymTracelessMultiplication : public MatrixBinaryOperator<R, T>
  {
  public:
    // Put public methods here. These should change very little over time.
    using MatrixBinaryOperator<R, T>::mR;
    using MatrixBinaryOperator<R, T>::mT;

    DEVICE_FUNCTION
    HermSymTracelessMultiplication(const R &pR, const T &pT) : MatrixBinaryOperator<R, T>(pR, pT) {}

    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<1> t1, Tag<1> t2) const
    {
      return getComponent(mR, 1_c, 1_c) * getComponent(mT, 1_c, 1_c) +
             getComponent(mR, 1_c, 2_c) * getComponent(mT, 2_c, 1_c) +
             getComponent(mR, 1_c, 3_c) * getComponent(mT, 3_c, 1_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<1> t1, Tag<2> t2) const
    {
      return getComponent(mR, 1_c, 1_c) * getComponent(mT, 1_c, 2_c) +
             getComponent(mR, 1_c, 2_c) * getComponent(mT, 2_c, 2_c) +
             getComponent(mR, 1_c, 3_c) * getComponent(mT, 3_c, 2_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<1> t1, Tag<3> t2) const
    {
      return getComponent(mR, 1_c, 1_c) * getComponent(mT, 1_c, 3_c) +
             getComponent(mR, 1_c, 2_c) * getComponent(mT, 2_c, 3_c) +
             getComponent(mR, 1_c, 3_c) * getComponent(mT, 3_c, 3_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<2> t1, Tag<1> t2) const
    {
      return getComponent(mR, 2_c, 1_c) * getComponent(mT, 1_c, 1_c) +
             getComponent(mR, 2_c, 2_c) * getComponent(mT, 2_c, 1_c) +
             getComponent(mR, 2_c, 3_c) * getComponent(mT, 3_c, 1_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<2> t1, Tag<2> t2) const
    {
      return getComponent(mR, 2_c, 1_c) * getComponent(mT, 1_c, 2_c) +
             getComponent(mR, 2_c, 2_c) * getComponent(mT, 2_c, 2_c) +
             getComponent(mR, 2_c, 3_c) * getComponent(mT, 3_c, 2_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<2> t1, Tag<3> t2) const
    {
      return getComponent(mR, 2_c, 1_c) * getComponent(mT, 1_c, 3_c) +
             getComponent(mR, 2_c, 2_c) * getComponent(mT, 2_c, 3_c) +
             getComponent(mR, 2_c, 3_c) * getComponent(mT, 3_c, 3_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<3> t1, Tag<1> t2) const
    {
      return getComponent(mR, 3_c, 1_c) * getComponent(mT, 1_c, 1_c) +
             getComponent(mR, 3_c, 2_c) * getComponent(mT, 2_c, 1_c) +
             getComponent(mR, 3_c, 3_c) * getComponent(mT, 3_c, 1_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<3> t1, Tag<2> t2) const
    {
      return getComponent(mR, 3_c, 1_c) * getComponent(mT, 1_c, 2_c) +
             getComponent(mR, 3_c, 2_c) * getComponent(mT, 2_c, 2_c) +
             getComponent(mR, 3_c, 3_c) * getComponent(mT, 3_c, 2_c);
    }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<3> t1, Tag<3> t2) const
    {
      return getComponent(mR, 3_c, 1_c) * getComponent(mT, 1_c, 3_c) +
             getComponent(mR, 3_c, 2_c) * getComponent(mT, 2_c, 3_c) +
             getComponent(mR, 3_c, 3_c) * getComponent(mT, 3_c, 3_c);
    }

    template <int N, int M> auto operator()(Tag<N> t1, Tag<M> t2) const
    {
      static_assert(N >= 1 && N <= 3 && M >= 1 && M <= 3, "Operator(): N and M must be between 1 and 3 for SymWrapper");
      return MatrixGet(t1, t2);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, std::decay_t<T> t, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
        DoEval::eval(t, idx...);
      }
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      auto cL = DoEval::eval(mR, idx...);
      auto cR = DoEval::eval(mT, idx...);
      device::array<decltype(cL[1] * cR[1]), 9> result;
      result[0] = cL[0] * cR[0] + cL[1] * cR[1] + cL[2] * cR[2];
      result[1] = cL[0] * cR[1] + cL[1] * cR[3] + cL[2] * cR[4];
      result[2] = cL[0] * cR[2] + cL[1] * cR[4] - cL[2] * (cR[0] + cR[3]);
      result[3] = conj(cL[1]) * cR[0] + cL[3] * cR[1] + cL[4] * cR[2];
      result[4] = conj(cL[1]) * cR[1] + cL[3] * cR[3] + cL[4] * cR[4];
      result[5] = conj(cL[1]) * cR[2] + cL[3] * cR[4] - cL[4] * (cR[0] + cR[3]);
      result[6] = conj(cL[2]) * cR[0] + conj(cL[4]) * cR[1] + cL[5] * cR[2];
      result[7] = conj(cL[2]) * cR[1] + conj(cL[4]) * cR[3] + cL[5] * cR[4];
      result[8] = conj(cL[2]) * cR[2] + conj(cL[4]) * cR[4] - cL[5] * (cR[0] + cR[3]);
      return result;
    }

    virtual std::string operatorString() const override { return "*"; }
  };

  template <typename R, typename T>
    requires(HasHermGet<R> && HasSymTracelessGet<T>)
  DEVICE_FORCEINLINE_FUNCTION auto operator*(const R &r, const T &t)
  {
    return HermSymTracelessMultiplication<R, T>(r, t);
  }

} // namespace TempLat

#endif
