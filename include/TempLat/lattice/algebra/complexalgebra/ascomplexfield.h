#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_ASCOMPLEXFIELD_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_ASCOMPLEXFIELD_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/operators/unaryoperator.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/complexfieldget.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/helpers/iscomplextype.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/parallel/device.h"

// IMPORTANT: this header must NOT include real.h / imag.h. They include this header (to add their
// third overload), so including them back would form a cycle. Everything used here (UnaryOperator,
// complexfieldget.h) is verified not to reach real.h / imag.h.

namespace TempLat
{
  /** @brief The real (N==0) or imaginary (N==1) component of a getter whose eval() returns a complex value.
   *
   * Derives from UnaryOperator so pre/postGet, confirmSpace, ghosts and toolBox forward to the child.
   **/
  template <typename R, int N> class ComplexEvalComponent : public UnaryOperator<R>
  {
    static_assert(N == 0 || N == 1, "ComplexEvalComponent: N must be 0 (real) or 1 (imaginary).");

  public:
    using UnaryOperator<R>::mR;

    ComplexEvalComponent(const R &pR) : UnaryOperator<R>(pR) {}

    static constexpr size_t NDim = GetNDim::get<R>();

    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      if constexpr (N == 0)
        return device::real(DoEval::eval(mR, idx...));
      else
        return device::imag(DoEval::eval(mR, idx...));
    }

    virtual std::string operatorString() const override { return N == 0 ? "Re " : "Im "; }
  };

  /** @brief Adapts a getter whose eval() returns a complex<C> value into a complex-field expression
   *  (i.e. one exposing ComplexFieldGet(Tag<0>)/ComplexFieldGet(Tag<1>)). This is the inverse of asFourier().
   *
   *  Deliberately a single node that evaluates its child ONCE and splits the result, rather than sugar for
   *  Complexify(Real(r), Imag(r)). PostGet::apply recurses structurally with no de-duplication, and a
   *  RandomGaussianField's `generation` counter is a refcount-shared host_ptr -- so two copies of the child
   *  in one tree would advance the RNG twice and draw the real and imaginary parts from different generations.
   *  Keep the child unique.
   *
   * Unit test: ctest -R test-ascomplexfield
   **/
  template <typename R> class AsComplexField : public UnaryOperator<R>
  {
  public:
    using UnaryOperator<R>::mR;

    AsComplexField(const R &pR) : UnaryOperator<R>(pR) {}

    static constexpr size_t NDim = GetNDim::get<R>();

    using ValueType = typename GetGetReturnType<std::decay_t<R>>::type; // complex<C>
    using ComponentType = typename ValueType::value_type;               // C

    auto ComplexFieldGet(Tag<0> t) const { return ComplexEvalComponent<R, 0>(mR); }
    auto ComplexFieldGet(Tag<1> t) const { return ComplexEvalComponent<R, 1>(mR); }

    template <int N> auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 1, "Operator(): N must be 0 or 1 for AsComplexField");
      return ComplexFieldGet(t);
    }

    /** @brief Returns device::array<C,2>, NOT complex<C>: the assignment loops index result[0]/result[1]. */
    template <typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto v = DoEval::eval(mR, idx...);
      device::array<ComponentType, 2> result;
      result[0] = device::real(v);
      result[1] = device::imag(v);
      return result;
    }

    std::string toString() const { return "asComplexField(" + GetString::get(mR) + ")"; }

    using Getter = ComplexFieldGetter;
    static constexpr size_t SHIFTIND = 0;
    static constexpr size_t size = 2;
  };

  /** @brief A getter that has an eval() returning a complex value, but does not itself speak the
   *  complex-field protocol. asComplexField() adapts exactly these. */
  template <typename T>
  concept HasComplexEval =
      HasEvalMethod<T> && !HasComplexFieldGet<T> && !IsComplexType<T> && GetGetReturnType<std::decay_t<T>>::isComplex;

  /**
   * @vocab-summary Reinterprets an expression whose eval() already returns a complex value as a complex field,
   * so the complex algebra applies to it.
   **/
  template <typename R>
    requires HasComplexEval<std::decay_t<R>>
  auto asComplexField(R &&r)
  {
    return AsComplexField<std::decay_t<R>>(std::forward<R>(r));
  }
} // namespace TempLat

#endif
