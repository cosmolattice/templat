#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMSYMTRACELESSMULTIPLYTRACE_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMSYMTRACELESSMULTIPLYTRACE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

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
  template <class R, class T> class SymSymTracelessMultiplicationTrace : public TempLat::BinaryOperator<R, T>
  {
  public:
    // Put public methods here. These should change very little over time.
    using TempLat::BinaryOperator<R, T>::mR;
    using TempLat::BinaryOperator<R, T>::mT;

    DEVICE_FUNCTION
    SymSymTracelessMultiplicationTrace(const R &pR, const T &pT) : TempLat::BinaryOperator<R, T>(pR, pT) {}

    DEVICE_FORCEINLINE_FUNCTION
    auto Get() const { return getComponent(mR, 1_c, 1_c) * getComponent(mT, 1_c, 1_c) +
                              getComponent(mR, 1_c, 2_c) * getComponent(mT, 2_c, 1_c) +
                              getComponent(mR, 1_c, 3_c) * getComponent(mT, 3_c, 1_c) +
                              getComponent(mR, 2_c, 1_c) * getComponent(mT, 1_c, 2_c) +
                              getComponent(mR, 2_c, 2_c) * getComponent(mT, 2_c, 2_c) +
                              getComponent(mR, 2_c, 3_c) * getComponent(mT, 3_c, 2_c) +
                              getComponent(mR, 3_c, 1_c) * getComponent(mT, 1_c, 3_c) +
                              getComponent(mR, 3_c, 2_c) * getComponent(mT, 2_c, 3_c) +
                              getComponent(mR, 3_c, 3_c) * getComponent(mT, 3_c, 3_c); }


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
        return cL[0] * cR[0] + cL[1] * cR[1] + cL[2] * cR[2] + cL[1] * cR[1] + cL[3] * cR[3] + cL[4] * cR[4] + cL[2] * cR[2] + cL[4] * cR[4] - cL[5] * (cR[0] + cR[3]);
    }

    virtual std::string operatorString() const override { return "*"; }
  };

  template <typename R, typename T>
    requires(HasSymGet<R> && HasSymTracelessGet<T> && !HasHermGet<R>)
  DEVICE_FORCEINLINE_FUNCTION auto multiplyTrace(const R &r, const T &t)
  {
    return SymSymTracelessMultiplicationTrace<R, T>(r, t);
  }

} // namespace TempLat

#endif
