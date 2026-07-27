#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2DOTTER_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2DOTTER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

#include "TempLat/lattice/algebra/operators/binaryoperator.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief Lie-algebra inner product <A,B> = sum_{a=1,2,3} A_a B_a of two SU(2)-valued expressions,
   *  returned as a scalar.
   *
   *  It evaluates each operand once through the FUSED eval path (one device::array<.,4> per operand)
   *  rather than contracting per component as sum_a A.SU2Get(a) * B.SU2Get(a). That distinction matters
   *  when an operand is a deep SU(2) expression - e.g. the clover/non-abelian magnetic field B4NA: its
   *  per-component SU2Get(a) expansion is the exponential decomposition of a product of plaquettes, and
   *  building it three times (a = 1,2,3) yields an expression object that overflows the host stack during
   *  construction. Evaluating each operand once keeps the object small and linear in the chain length.
   *
   *  Being scalar-valued, it plugs straight into reductions such as average(). Component 0 (the SU(2)
   *  identity/trace part) is intentionally dropped, matching the Lie-algebra convention of
   *  SU2LieAlgebraField.
   *
   * Unit test: ctest -R test-su2dotter
   **/
  template <typename R, typename T> class SU2Dotter : public BinaryOperator<R, T>
  {
  public:
    using BinaryOperator<R, T>::mR;
    using BinaryOperator<R, T>::mT;

    SU2Dotter(const R &pR, const T &pT) : BinaryOperator<R, T>(pR, pT) {}

    template <typename... IDX>
      requires requires(std::decay_t<R> r, std::decay_t<T> t, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
        DoEval::eval(t, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto rr = DoEval::eval(mR, idx...);
      const auto tt = DoEval::eval(mT, idx...);
      return rr[1] * tt[1] + rr[2] * tt[2] + rr[3] * tt[3];
    }

    virtual std::string operatorString() const override { return "·"; }
  };

  /**
   * @vocab-summary Lie-algebra inner product $\langle A,B\rangle = \sum_{a=1}^{3} A_a B_a$, evaluated in one
   * fused pass rather than by building three component expressions.
   **/
  template <class R, class T>
    requires(HasSU2Get<R> && HasSU2Get<T>)
  auto su2dotter(const R &r, const T &t)
  {
    return SU2Dotter<R, T>(r, t);
  }
} // namespace TempLat

#endif
