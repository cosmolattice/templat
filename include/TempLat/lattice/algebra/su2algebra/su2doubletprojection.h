#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2DOUBLETPROJECTION_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_SU2DOUBLETPROJECTION_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/complexalgebra/real.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/su2algebra/su2binaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletdotter.h"
#include "TempLat/lattice/algebra/su2algebra/su2generators.h"
#include "TempLat/lattice/algebra/su2algebra/su2su2doubletmultiply.h"
#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief The SU(2) element with components Re<x, i sigma_b c>, b = 0 ... 3, for doublets x and c.
   *
   *  Equals MakeSU2(b, Real(scalar_prod(x, i_sigma(b) * c))), which evaluates x and c once per component; this
   *  evaluates them once. Matters when x is costly, e.g. a covariant shift U(x) phi(x + i) in a matter current.
   *
   *  Unit test: ctest -R test-su2doubletprojection
   **/
  template <typename R, typename T> class SU2DoubletProjection : public SU2BinaryOperator<R, T>
  {
  public:
    using SU2BinaryOperator<R, T>::mR;
    using SU2BinaryOperator<R, T>::mT;

    SU2DoubletProjection(const R &pR, const T &pT) : SU2BinaryOperator<R, T>(pR, pT) {}

    template <int B> auto SU2Get(Tag<B> b) const { return Real(scalar_prod(mR, Constants::i_sigma(b) * mT)); }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, std::decay_t<T> t, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
        DoEval::eval(t, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto x = DoEval::eval(mR, idx...);
      const auto c = DoEval::eval(mT, idx...);
      // i sigma_b c as in SU2SU2DoubletMultiplication with the unit element b, dotted with x as in SU2DoubletDotter.
      device::array<std::decay_t<decltype(x[0] * c[0])>, 4> result;
      result[0] = x[0] * c[0] + x[1] * c[1] + x[2] * c[2] + x[3] * c[3];
      result[1] = -x[0] * c[3] + x[1] * c[2] - x[2] * c[1] + x[3] * c[0];
      result[2] = x[0] * c[2] + x[1] * c[3] - x[2] * c[0] - x[3] * c[1];
      result[3] = -x[0] * c[1] + x[1] * c[0] + x[2] * c[3] - x[3] * c[2];
      return result;
    }

    virtual std::string operatorString() const override { return "·iσ"; }
  };

  /**
   * @vocab-summary The SU(2) element $\mathrm{Re}\langle x, i\sigma_b c\rangle$, $b = 0 \ldots 3$, of two doublets,
   * e.g. the bilinear in the SU(2) current. Evaluates $x$ and $c$ once.
   * @vocab-signature doubletProjection(x, c)
   **/
  template <typename R, typename T> SU2DoubletProjection<R, T> doubletProjection(const R &x, const T &c)
  {
    return {x, c};
  }
} // namespace TempLat

#endif
