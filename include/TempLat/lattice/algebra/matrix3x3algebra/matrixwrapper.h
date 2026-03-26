#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_MATRIXWRAPPER_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_MATRIXWRAPPER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/allmatrixoperator.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"

namespace TempLat
{
  /** @brief A class which wraps two objects as a Matrix3x3  field.
   *
   *
   * Unit test: ctest -R test-complexwrapper
   **/
  template <typename R11, typename R12, typename R13, typename R21, typename R22, typename R23, typename R31, typename R32, typename R33> class MatrixWrapper : public MatrixOperator
  {
  public:
    // Put public methods here. These should change very little over time.

    DEVICE_FUNCTION
    MatrixWrapper() = default;

    DEVICE_FUNCTION
    MatrixWrapper(const R11 &pR11, const R12 &pR12,const R13 &pR13, const R21 &pR21, const R22 &pR22, const R23 &pR23, const R31 &pR31, const R32 &pR32, const R33 &pR33):
    mR11(pR11),
    mR12(pR12),
    mR13(pR13),
    mR21(pR21),
    mR22(pR22),
    mR23(pR23),
    mR31(pR31),
    mR32(pR32),
    mR33(pR33)
    {}

    DEVICE_FUNCTION
    MatrixWrapper(const MatrixWrapper &) = default;

    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<0> t) const { return mR11; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<1> t) const { return mR12; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<2> t) const { return mR13; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<3> t) const { return mR21; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<4> t) const { return mR22; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<5> t) const { return mR23; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<6> t) const { return mR31; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<7> t) const { return mR32; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<8> t) const { return mR33; }

    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<1> t1, Tag<1> t2) const { return mR11; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<1> t1, Tag<2> t2) const { return mR12; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<1> t1, Tag<3> t2) const { return mR13; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<2> t1, Tag<1> t2) const { return mR21; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<2> t1, Tag<2> t2) const { return mR22; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<2> t1, Tag<3> t2) const { return mR23; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<3> t1, Tag<1> t2) const { return mR31; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<3> t1, Tag<2> t2) const { return mR32; }
    DEVICE_FORCEINLINE_FUNCTION
    auto MatrixGet(Tag<3> t1, Tag<3> t2) const { return mR33; }

    template <int N> DEVICE_FORCEINLINE_FUNCTION auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 8, "Operator(): N must be between 0 and 8 for MatrixWrapper");
      return MatrixGet(t);
    }

    template <int N, int M> DEVICE_FORCEINLINE_FUNCTION auto operator()(Tag<N> t1, Tag<M> t2) const
    {
      static_assert(N >= 1 && N <= 3 && M >= 1 && M <= 3, "Operator(): N and M must be between 1 and 3 for MatrixWrapper");
      return MatrixGet(t1, t2);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R11> r11, std::decay_t<R12> r12, std::decay_t<R13> r13, std::decay_t<R21> r21, std::decay_t<R22> r22, std::decay_t<R23> r23, std::decay_t<R31> r31, std::decay_t<R32> r32, std::decay_t<R33> r33, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r11, idx...);
        DoEval::eval(r12, idx...);
        DoEval::eval(r13, idx...);
        DoEval::eval(r21, idx...);
        DoEval::eval(r22, idx...);
        DoEval::eval(r23, idx...);
        DoEval::eval(r31, idx...);
        DoEval::eval(r32, idx...);
        DoEval::eval(r33, idx...);
      }
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<decltype(DoEval::eval(mR11, idx...)), 9> result;
      result[0] = DoEval::eval(mR11, idx...);
      result[1] = DoEval::eval(mR12, idx...);
      result[2] = DoEval::eval(mR13, idx...);
      result[3] = DoEval::eval(mR21, idx...);
      result[4] = DoEval::eval(mR22, idx...);
      result[5] = DoEval::eval(mR23, idx...);
      result[6] = DoEval::eval(mR31, idx...);
      result[7] = DoEval::eval(mR32, idx...);
      result[8] = DoEval::eval(mR33, idx...);
      return result;
    }

    void preGet()
    {
      PreGet::apply(mR11);
      PreGet::apply(mR12);
      PreGet::apply(mR13);
      PreGet::apply(mR21);
      PreGet::apply(mR22);
      PreGet::apply(mR23);
      PreGet::apply(mR31);
      PreGet::apply(mR32);
      PreGet::apply(mR33);
    }

    void postGet()
    {
      PostGet::apply(mR11);
      PostGet::apply(mR12);
      PostGet::apply(mR13);
      PostGet::apply(mR21);
      PostGet::apply(mR22);
      PostGet::apply(mR23);
      PostGet::apply(mR31);
      PostGet::apply(mR32);
      PostGet::apply(mR33);
    }

    std::string toString() const { return "Matrix( (" + GetString::get(mR11) + "," + GetString::get(mR12) + "," + GetString::get(mR13) + ") , (" + GetString::get(mR21) + "," + GetString::get(mR22) + "," + GetString::get(mR23) + ") , (" + GetString::get(mR31) + "," + GetString::get(mR32) + "," + GetString::get(mR33) + ") "; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    R11 mR11;
    R12 mR12;
    R13 mR13;
    R21 mR21;
    R22 mR22;
    R23 mR23;
    R31 mR31;
    R32 mR32;
    R33 mR33;
  };

  template <typename R11, typename R12, typename R13, typename R21, typename R22, typename R23, typename R31, typename R32, typename R33>
  DEVICE_FORCEINLINE_FUNCTION MatrixWrapper<R11, R12, R13, R21, R22, R23, R31, R32, R33> ConstructMatrix3x3(const R11 &r11, const R12 &r12, const R13 &r13, const R21 &r21, const R22 &r22, const R23 &r23, const R31 &r31, const R32 &r32, const R33 &r33)
  {
    return {r11, r12, r13, r21, r22, r23, r31, r32, r33};
  }
} // namespace TempLat

#endif
