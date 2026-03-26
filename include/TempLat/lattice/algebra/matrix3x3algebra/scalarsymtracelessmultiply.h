#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SCALARSYMTRACELESSMULTIPLY_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SCALARSYMTRACELESSMULTIPLY_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/helpers/hasstaticgetter.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessbinaryoperator.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/isscalartype.h"
#include <type_traits>

namespace TempLat
{
  /** @brief A class which implements scalar multiplication over symmetric traceless numbers.
   *
   * Unit test: ctest -R test-scalarcomplexfieldmultiply
   **/

  template <typename R, typename T> class ScalarSymTracelessMultiply : public SymTracelessBinaryOperator<R, T>
  {
  public:
    // Put public methods here. These should change very little over time.

    using SymTracelessBinaryOperator<R, T>::mR;
    using SymTracelessBinaryOperator<R, T>::mT;

    DEVICE_FUNCTION
    ScalarSymTracelessMultiply(const R &pR, const T &pT) : SymTracelessBinaryOperator<R, T>(pR, pT) {}

    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<0> t) const { return mR * getComponent(mT, 0_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t) const { return mR * getComponent(mT, 1_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t) const { return mR * getComponent(mT, 2_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t) const { return mR * getComponent(mT, 3_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<4> t) const { return mR * getComponent(mT, 4_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<5> t) const { return mR * (- getComponent(mT, 0_c) - getComponent(mT, 3_c)); }

    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t1, Tag<1> t2) const { return SymTracelessGet(0_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t1, Tag<2> t2) const { return SymTracelessGet(1_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t1, Tag<3> t2) const { return SymTracelessGet(2_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t1, Tag<1> t2) const { return SymTracelessGet(1_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t1, Tag<2> t2) const { return SymTracelessGet(3_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t1, Tag<3> t2) const { return SymTracelessGet(4_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t1, Tag<1> t2) const { return SymTracelessGet(2_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t1, Tag<2> t2) const { return SymTracelessGet(4_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t1, Tag<3> t2) const { return SymTracelessGet(5_c); }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, std::decay_t<T> t, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
        DoEval::eval(t, idx...);
      }
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto symtraceless = DoEval::eval(mT, idx...);
      const auto scalar = DoEval::eval(mR, idx...);
      device::array<decltype(scalar * symtraceless[0]), 5> result;
      result[0] = scalar * symtraceless[0];
      result[1] = scalar * symtraceless[1];
      result[2] = scalar * symtraceless[2];
      result[3] = scalar * symtraceless[3];
      result[4] = scalar * symtraceless[4];
      return result;
    }

    virtual std::string operatorString() const override { return "*"; }
  };

  template <typename R, typename T>
    requires(IsScalarType<R> && HasSymTracelessGet<T>)
  DEVICE_FORCEINLINE_FUNCTION auto operator*(const R &r, const T &t)
  {
    return ScalarSymTracelessMultiply<R, T>(r, t);
  }

  template <typename R, typename T>
    requires(HasSymTracelessGet<R> && IsScalarType<T>)
  DEVICE_FORCEINLINE_FUNCTION auto operator*(const R &r, const T &t)
  {
    return ScalarSymTracelessMultiply<T, R>(t, r);
  }

  template <typename R, typename T>
    requires(HasSymTracelessGet<R> && IsScalarType<T>)
  DEVICE_FORCEINLINE_FUNCTION auto operator/(const R &r, const T &t)
  {
    return ScalarSymTracelessMultiply(1_c / t, r);
  }
} // namespace TempLat

#endif
