#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMSYMTRACELESSMULTIPLY_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMSYMTRACELESSMULTIPLY_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/matrixbinaryoperator.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/symget.h"
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
  template <class R, class T> class SymSymTracelessMultiplication : public MatrixBinaryOperator<R, T>
  {
  public:
    // Put public methods here. These should change very little over time.
    using MatrixBinaryOperator<R, T>::mR;
    using MatrixBinaryOperator<R, T>::mT;

    SymSymTracelessMultiplication(const R &pR, const T &pT) : MatrixBinaryOperator<R, T>(pR, pT) {}

    auto MatrixGet(Tag<1> t1, Tag<1> t2) const
    {
      return getComponent(mR, 1_c, 1_c) * getComponent(mT, 1_c, 1_c) +
             getComponent(mR, 1_c, 2_c) * getComponent(mT, 2_c, 1_c) +
             getComponent(mR, 1_c, 3_c) * getComponent(mT, 3_c, 1_c);
    }
    auto MatrixGet(Tag<1> t1, Tag<2> t2) const
    {
      return getComponent(mR, 1_c, 1_c) * getComponent(mT, 1_c, 2_c) +
             getComponent(mR, 1_c, 2_c) * getComponent(mT, 2_c, 2_c) +
             getComponent(mR, 1_c, 3_c) * getComponent(mT, 3_c, 2_c);
    }
    auto MatrixGet(Tag<1> t1, Tag<3> t2) const
    {
      return getComponent(mR, 1_c, 1_c) * getComponent(mT, 1_c, 3_c) +
             getComponent(mR, 1_c, 2_c) * getComponent(mT, 2_c, 3_c) +
             getComponent(mR, 1_c, 3_c) * getComponent(mT, 3_c, 3_c);
    }
    auto MatrixGet(Tag<2> t1, Tag<1> t2) const
    {
      return getComponent(mR, 2_c, 1_c) * getComponent(mT, 1_c, 1_c) +
             getComponent(mR, 2_c, 2_c) * getComponent(mT, 2_c, 1_c) +
             getComponent(mR, 2_c, 3_c) * getComponent(mT, 3_c, 1_c);
    }
    auto MatrixGet(Tag<2> t1, Tag<2> t2) const
    {
      return getComponent(mR, 2_c, 1_c) * getComponent(mT, 1_c, 2_c) +
             getComponent(mR, 2_c, 2_c) * getComponent(mT, 2_c, 2_c) +
             getComponent(mR, 2_c, 3_c) * getComponent(mT, 3_c, 2_c);
    }
    auto MatrixGet(Tag<2> t1, Tag<3> t2) const
    {
      return getComponent(mR, 2_c, 1_c) * getComponent(mT, 1_c, 3_c) +
             getComponent(mR, 2_c, 2_c) * getComponent(mT, 2_c, 3_c) +
             getComponent(mR, 2_c, 3_c) * getComponent(mT, 3_c, 3_c);
    }
    auto MatrixGet(Tag<3> t1, Tag<1> t2) const
    {
      return getComponent(mR, 3_c, 1_c) * getComponent(mT, 1_c, 1_c) +
             getComponent(mR, 3_c, 2_c) * getComponent(mT, 2_c, 1_c) +
             getComponent(mR, 3_c, 3_c) * getComponent(mT, 3_c, 1_c);
    }
    auto MatrixGet(Tag<3> t1, Tag<2> t2) const
    {
      return getComponent(mR, 3_c, 1_c) * getComponent(mT, 1_c, 2_c) +
             getComponent(mR, 3_c, 2_c) * getComponent(mT, 2_c, 2_c) +
             getComponent(mR, 3_c, 3_c) * getComponent(mT, 3_c, 2_c);
    }
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
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      auto cL = DoEval::eval(mR, idx...);
      auto cR = DoEval::eval(mT, idx...);
      device::array<decltype(cL[0] * cR[0]), 9> result;
      result[0] = cL[0] * cR[0] + cL[1] * cR[1] + cL[2] * cR[2];
      result[1] = cL[0] * cR[1] + cL[1] * cR[3] + cL[2] * cR[4];
      result[2] = cL[0] * cR[2] + cL[1] * cR[4] - cL[2] * (cR[0] + cR[3]);
      result[3] = cL[1] * cR[0] + cL[3] * cR[1] + cL[4] * cR[2];
      result[4] = cL[1] * cR[1] + cL[3] * cR[3] + cL[4] * cR[4];
      result[5] = cL[1] * cR[2] + cL[3] * cR[4] - cL[4] * (cR[0] + cR[3]);
      result[6] = cL[2] * cR[0] + cL[4] * cR[1] + cL[5] * cR[2];
      result[7] = cL[2] * cR[1] + cL[4] * cR[3] + cL[5] * cR[4];
      result[8] = cL[2] * cR[2] + cL[4] * cR[4] - cL[5] * (cR[0] + cR[3]);
      return result;
    }

    virtual std::string operatorString() const override { return "*"; }
  };

  template <typename R, typename T>
    requires(HasSymGet<R> && HasSymTracelessGet<T>)
  auto operator*(const R &r, const T &t)
  {
    return SymSymTracelessMultiplication<R, T>(r, t);
  }

} // namespace TempLat

#endif
