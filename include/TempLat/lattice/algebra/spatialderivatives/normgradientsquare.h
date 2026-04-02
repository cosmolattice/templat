#ifndef TEMPLAT_LATTICE_ALGEBRA_SPATIALDERIVATIVES_NORMGRADIENTSQUARE_H
#define TEMPLAT_LATTICE_ALGEBRA_SPATIALDERIVATIVES_NORMGRADIENTSQUARE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler  Year: 2025

#include "TempLat/lattice/algebra/spatialderivatives/latticeforwardgradient.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/algebra/operators/unaryoperator.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/util/constexpr_for.h"
#include "TempLat/util/tuple_tools.h"

#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A class which computes the gradient square norm.
   *
   *
   * Unit test: ctest -R test-normgradientsquare
   **/
  template <typename R> class NormGradientSquare : public UnaryOperator<R>
  {
  public:
    // Put public methods here. These should change very little over time.
    using GetReturnType = typename GetGetReturnType<R>::type;
    using FloatType = typename GetFloatType<GetReturnType>::type;
    static constexpr size_t NDim = GetNDim::get<R>();

    using UnaryOperator<R>::mR;

    DEVICE_FUNCTION
    NormGradientSquare(const R &pR) : UnaryOperator<R>(pR), dx2(pow<2>(GetDx::getDx(pR))) {}

    std::string toString() const { return "|Grad(" + GetString::get(mR) + ")|^2"; }

    void doWeNeedGhosts() const { mR.confirmGhostsUpToDate(); }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      if constexpr (UnaryOperator<R>::getNDim() == 0)
        return ZeroType();
      else {
        const auto midval = DoEval::eval(mR, idx...);
        FloatType result{};
        constexpr_for<0, NDim>([&](const auto _d) {
          constexpr size_t d = decltype(_d)::value;
          device::apply(
              [&](const auto &...shifted_idx) { result += pow<2>(DoEval::eval(mR, shifted_idx...) - midval); },
              tuple_add_to_nth<d, 1>(device::tie(idx...)));
        });
        return result / dx2;
      }
    }

    template <typename S> DEVICE_INLINE_FUNCTION auto d(const S &other)
    {
      return 2 * LatForwardGrad(mR) * LatForwardGrad(mR.d(other));
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    FloatType dx2;
  };

  template <size_t NDim_ = 0, typename R>
    requires HasEvalMethod<R>
  DEVICE_INLINE_FUNCTION auto Grad2(R pR)
  {
    static_assert(NDim_ == 0 || NDim_ == GetNDim::get<R>(),
                  "Explicit NDim does not match the NDim deduced from expression type R.");
    return NormGradientSquare<R>(pR);
  }

  template <size_t NDim_ = 0, typename R>
    requires(!HasEvalMethod<R>)
  DEVICE_INLINE_FUNCTION auto Grad2(R pR)
  {
    return ZeroType();
  }
} // namespace TempLat

#endif
