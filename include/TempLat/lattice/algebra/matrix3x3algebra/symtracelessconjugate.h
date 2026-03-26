#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSCONJUGATE_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSCONJUGATE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessunaryoperator.h"
#include "TempLat/lattice/algebra/helpers/hasstaticgetter.h"
#include "TempLat/lattice/algebra/operators/complexconjugate.h"

namespace TempLat
{

  /** @brief A class which implements conjugation of all the elements of a symmetric-traceless matrix.
   *
   * Unit test: ctest -R test-symtracelessconjugate
   **/
  template <typename R> class SymTracelessConjugate : public SymTracelessUnaryOperator<R>
  {
  public:
    // Put public methods here. These should change very little over time.

    using SymTracelessUnaryOperator<R>::mR;

    DEVICE_FUNCTION
    SymTracelessConjugate(const R &pR) : SymTracelessUnaryOperator<R>(pR) {}

    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<0> t) const { return conj(getComponent(mR, 0_c)); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t) const { return conj(getComponent(mR, 1_c)); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t) const { return conj(getComponent(mR, 2_c)); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t) const { return conj(getComponent(mR, 3_c)); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<4> t) const { return conj(getComponent(mR, 4_c)); }

    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t1, Tag<1> t2) const { return (*this).SymTracelessGet(0_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t1, Tag<2> t2) const { return (*this).SymTracelessGet(1_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<1> t1, Tag<3> t2) const { return (*this).SymTracelessGet(2_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t1, Tag<1> t2) const { return (*this).SymTracelessGet(1_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t1, Tag<2> t2) const { return (*this).SymTracelessGet(3_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<2> t1, Tag<3> t2) const { return (*this).SymTracelessGet(4_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t1, Tag<1> t2) const { return (*this).SymTracelessGet(2_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t1, Tag<2> t2) const { return (*this).SymTracelessGet(4_c); }
    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t1, Tag<3> t2) const { return - (*this).SymTracelessGet(0_c) - (*this).SymTracelessGet(3_c); }

    template <int N, int M> DEVICE_FORCEINLINE_FUNCTION auto operator()(Tag<N> t1, Tag<M> t2) const
    {
      static_assert(N >= 1 && N <= 3 && M >= 1 && M <= 3, "Operator(): N and M must be between 1 and 3 for SymWrapper");
      return SymTracelessGet(t1, t2);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      auto cL = DoEval::eval(mR, idx...);
      device::array<decltype(cL[0]), 5> result; //TODE: Jorge: For now this assummes cL[0] is already a complex. I only plan to use this in Fourier space, where i know it is, but this should no taken for granted.
      result[0] = conj(cL[0]);
      result[1] = conj(cL[1]);
      result[2] = conj(cL[2]);
      result[3] = conj(cL[3]);
      result[4] = conj(cL[4]);
      return result;
    }

    virtual std::string operatorString() const override { return "+"; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };

  template <typename R>
  requires HasSymTracelessGet<R>
  DEVICE_FORCEINLINE_FUNCTION auto conj(const R &r)
  {
    return SymTracelessConjugate<R>(r);
  }

} // namespace TempLat

#endif
