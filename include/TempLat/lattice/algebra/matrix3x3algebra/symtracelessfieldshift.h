#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSFIELDSHIFT_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSFIELDSHIFT_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymtracelessget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessunaryoperator.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

namespace TempLat
{
  /** @brief A class which implements spatial shifts for symmetric traceless algebra.
   *
   *
   * Unit test: ctest -R test-complexfieldshift
   **/
  template <typename R, int... N> class SymTracelessFieldShifter : public SymTracelessUnaryOperator<R>
  {
  public:
    // Put public methods here. These should change very little over time.
    static constexpr size_t dim = sizeof...(N);
    static constexpr auto shifts = device::make_tuple(N...);

    using SymTracelessUnaryOperator<R>::mR;

    SymTracelessFieldShifter(const R &pR) : SymTracelessUnaryOperator<R>(pR)
    {
      shiftString = shift<N...>(mR.SymTracelessFieldShifter(0_c)).getString({N...});
    }

    template <int M> DEVICE_FORCEINLINE_FUNCTION auto SymTracelessGet(Tag<M> t) const
    {
      static_assert(M >= 0 && M <= 4, "SymTracelessFieldShifter: M must be 0-4 for SymTracelessFieldShifter");
      return shift<N...>(mR.SymTracelessGet(t));
    }

    template <int M1, int M2> DEVICE_FORCEINLINE_FUNCTION auto SymTracelessGet(Tag<M1> t1, Tag<M2> t2) const
    {
      static_assert(M1 >= 1 && M1 <= 3 && M2 >= 1 && M2 <= 3, "SymTracelessFieldShifter: M1 and M2 must be 0-2 for SymTracelessFieldShifter");
      return shift<N...>(mR.SymTracelessGet(t1, t2));
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      auto tup = device::tie(idx...);
      constexpr_for<0, dim>([&](const auto _d) {
        constexpr size_t d = decltype(_d)::value;
        tup = tuple_add_to_nth<d, device::get<d>(shifts)>(tup);
      });
      return device::apply([&](const auto &...shifted_idx) { return DoEval::eval(mR, shifted_idx...); }, tup);
    }

    std::string operatorString() const { return shiftString; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    std::string shiftString;
  };
  template <typename R, int _N> class SymTracelessFieldShifterByOne : public SymTracelessUnaryOperator<R>
  {
    static_assert(_N != 0, "_N cannot be 0.");

    static constexpr int N = _N > 0 ? _N : -_N;
    static constexpr int dir = _N > 0 ? 1 : -1;

  public:
    // Put public methods here. These should change very little over time.

    using SymTracelessUnaryOperator<R>::mR;

    SymTracelessFieldShifterByOne(const R &pR) : SymTracelessUnaryOperator<R>(pR) {}

    template <int M> DEVICE_FORCEINLINE_FUNCTION auto SymTracelessGet(Tag<M> t) const
    {
      static_assert(M >= 0 && M <= 4, "SymTracelessFieldShifterByOne: M must be 0-4 for SymTracelessGet");
      return shift<_N>(mR.SymTracelessGet(t));
    }

    template <int M1, int M2> DEVICE_FORCEINLINE_FUNCTION auto SymTracelessGet(Tag<M1> t1, Tag<M2> t2) const
    {
      static_assert(M1 >= 1 && M1 <= 3 && M2 >= 1 && M2 <= 3, "SymTracelessFieldShifterByOne: M1 and M2 must be 0-2 for SymTracelessGet");
      return shift<_N>(mR.SymTracelessGet(t1, t2));
    }

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return device::apply([&](const auto &...shifted_idx) { return DoEval::eval(mR, shifted_idx...); },
                           tuple_add_to_nth<N - 1, dir>(device::tie(idx...)));
    }

    std::string toString() const { return GetString::get(mR) + "(->" + std::to_string(N) + ")"; }
  };

  template <int... shifts, class R>
    requires((sizeof...(shifts) > 1) && HasSymTracelessGet<R>)
  auto shift(const R &pR)
  {
    return SymTracelessFieldShifter<R, shifts...>(pR);
  }

  template <int N, class R>
    requires HasSymTracelessGet<R>
  auto shift(const R &pR)
  {
    return SymTracelessFieldShifterByOne<R, N>(pR);
  }

  template <class R, int N>
    requires HasSymTracelessGet<R>
  auto shift(const R &pR, Tag<N> t)
  {
    return SymTracelessFieldShifterByOne<R, N>(pR);
  }
} // namespace TempLat

#endif
