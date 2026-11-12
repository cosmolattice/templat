#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNMENT_WHOLELIST_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNMENT_WHOLELIST_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <type_traits>
#include <utility>

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/istemplatgettable.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/listoperators/listadd.h"
#include "TempLat/lattice/algebra/listoperators/listmultiply.h"
#include "TempLat/lattice/algebra/listoperators/listsubtract.h"
#include "TempLat/lattice/algebra/listoperators/listunaryminus.h"
#include "TempLat/lattice/field/assignment/assignviews.h"
#include "TempLat/parallel/device.h"
#include "TempLat/util/constexpr_for.h"

namespace TempLat
{
  /** @brief Evaluating a list right-hand side for all members at once ("whole"), so that members can share work, see
   *  plaqGrad. wholeOf(list) converts a list expression on the host: lists offering wholeList() share; sums,
   *  differences, negations and real multiples of lists combine member values component by component (widening e.g.
   *  integer zeros or float to double); any other list is evaluated member by member. A whole evaluation has size,
   *  evalWhole(idx...) (the array of member values) and shares (whether any part of it shares work).
   *
   *  Unit test: ctest -R test-plaquettegradient
   **/
  namespace WholeListDetail
  {
    /** @brief Component-wise a + s * b (s = +-1), for scalar or array values; b's type widens to a's (WidensTo). */
    template <typename V, typename W> DEVICE_INLINE_FUNCTION V combine(const V &a, const W &b, int s)
    {
      if constexpr (std::is_arithmetic_v<V>)
        return a + s * b;
      else {
        V out;
        for (size_t c = 0; c < a.size(); ++c)
          out[c] = a[c] + s * b[c];
        return out;
      }
    }
    template <typename V, typename S> DEVICE_INLINE_FUNCTION V scale(const S &s, const V &a)
    {
      if constexpr (std::is_arithmetic_v<V>)
        return s * a;
      else {
        V out;
        for (size_t c = 0; c < a.size(); ++c)
          out[c] = s * a[c];
        return out;
      }
    }

    /** @brief The members of any list, evaluated one by one. */
    template <typename T, typename... Ms> struct Members {
      static constexpr size_t size = sizeof...(Ms);
      static constexpr bool shares = false;
      device::tuple<Ms...> members;

      template <typename... IDX> DEVICE_INLINE_FUNCTION auto evalWhole(const IDX &...idx) const
      {
        return device::apply(
            [&](const auto &...m) {
              using V = std::common_type_t<std::remove_cvref_t<decltype(DoEval::eval(m, idx...))>...>;
              using E = std::conditional_t<std::is_arithmetic_v<V>, T, V>;
              return device::array<E, sizeof...(Ms)>{E(DoEval::eval(m, idx...))...};
            },
            members);
      }
      void preGet()
      {
        device::apply([](auto &...m) { (PreGet::apply(m), ...); }, members);
      }
      void postGet()
      {
        device::apply([](auto &...m) { (PostGet::apply(m), ...); }, members);
      }
    };

    /** @brief A + s * B member by member. */
    template <typename A, typename B> struct Sum {
      static constexpr size_t size = A::size;
      static constexpr bool shares = A::shares || B::shares;
      A a;
      B b;
      int s;

      template <typename... IDX> DEVICE_INLINE_FUNCTION auto evalWhole(const IDX &...idx) const
      {
        auto out = a.evalWhole(idx...);
        const auto other = b.evalWhole(idx...);
        for (size_t m = 0; m < size; ++m)
          out[m] = combine(out[m], other[m], s);
        return out;
      }
      void preGet()
      {
        a.preGet();
        b.preGet();
      }
      void postGet()
      {
        a.postGet();
        b.postGet();
      }
    };

    /** @brief The real scalar expression S times every member of A. */
    template <typename S, typename A> struct Scaled {
      static constexpr size_t size = A::size;
      static constexpr bool shares = A::shares;
      S scalar;
      A a;

      template <typename... IDX> DEVICE_INLINE_FUNCTION auto evalWhole(const IDX &...idx) const
      {
        const auto s = DoEval::eval(scalar, idx...);
        auto out = a.evalWhole(idx...);
        for (size_t m = 0; m < size; ++m)
          out[m] = scale(s, out[m]);
        return out;
      }
      void preGet()
      {
        PreGet::apply(scalar);
        a.preGet();
      }
      void postGet()
      {
        PostGet::apply(scalar);
        a.postGet();
      }
    };

    template <typename L>
    concept IsList = IsTempLatGettable<0, L>;

    template <typename V> struct ElementOf {
      using type = V;
    };
    template <typename E, size_t K> struct ElementOf<device::array<E, K>> {
      using type = E;
      static constexpr size_t extent = K;
    };
    template <typename V> constexpr size_t extentOf()
    {
      if constexpr (requires { ElementOf<V>::extent; })
        return ElementOf<V>::extent;
      else
        return 0;
    }

    /** @brief Member values of type From add into To without loss: both scalars, or arrays of equal size, with
     *  arithmetic elements whose common type is To's. E.g. the integer zeros of MakeSU2(a, ZeroType()) into SU(2)
     *  values. */
    template <typename From, typename To>
    concept WidensTo = std::is_arithmetic_v<typename ElementOf<From>::type> &&
                       std::is_arithmetic_v<typename ElementOf<To>::type> && extentOf<From>() == extentOf<To>() &&
                       std::is_same_v<std::common_type_t<typename ElementOf<From>::type, typename ElementOf<To>::type>,
                                      typename ElementOf<To>::type>;

    template <typename T, typename L, size_t... m> auto members(L &l, std::index_sequence<m...>)
    {
      using Getter = typename std::remove_cvref_t<L>::Getter;
      return Members<T, std::remove_cvref_t<decltype(Getter::get(l, Tag<int(m)>()))>...>{
          device::make_tuple(Getter::get(l, Tag<int(m)>())...)};
    }

    using AssignViewsDetail::AnIdx;
    template <typename W, size_t... d>
    auto memberValue(std::index_sequence<d...>)
        -> std::remove_cvref_t<decltype(std::declval<const W &>().evalWhole(std::declval<AnIdx<d>>()...)[0])>;
    template <typename S, size_t... d>
    auto scalarValue(std::index_sequence<d...>)
        -> std::remove_cvref_t<decltype(DoEval::eval(std::declval<const S &>(), std::declval<AnIdx<d>>()...))>;

    /** @brief The member value type of a whole evaluation, and the value type of a scalar expression, at an
     *  NDim-dimensional site. */
    template <typename W, size_t NDim> using MemberValue = decltype(memberValue<W>(std::make_index_sequence<NDim>()));
    template <typename S, size_t NDim> using ScalarValue = decltype(scalarValue<S>(std::make_index_sequence<NDim>()));
  } // namespace WholeListDetail

  /** @brief Whether the whole evaluation of a list of type L shares work, i.e. contains a list offering wholeList.
   *  Decided from the type alone, so that lists without one never build a whole evaluation. */
  template <typename L> constexpr bool wholeListShares()
  {
    using nakedL = std::remove_cvref_t<L>;
    if constexpr (requires(const nakedL &l) { l.wholeList(); })
      return true;
    else if constexpr (
        requires { typename nakedL::IsWholeListSum; } || requires { typename nakedL::IsWholeListProduct; })
      return wholeListShares<decltype(std::declval<nakedL &>().mR)>() ||
             wholeListShares<decltype(std::declval<nakedL &>().mT)>();
    else if constexpr (requires { typename nakedL::IsWholeListNegation; })
      return wholeListShares<decltype(std::declval<nakedL &>().mR)>();
    else
      return false;
  }

  /** @brief A whole evaluation as the right-hand side of a collection assignment: its value at a site is the array
   *  of all member values, see AssignKernel. */
  template <typename W> struct WholeAssignment {
    static constexpr bool fusesComponents = true;
    static constexpr size_t size = W::size;
    W whole;

    template <typename... IDX> DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return whole.evalWhole(idx...);
    }
    void preGet() { whole.preGet(); }
    void postGet() { whole.postGet(); }
  };

  /** @brief The whole evaluation of the list l of size N with NDim-dimensional members; T is the field type. */
  template <typename T, size_t NDim, size_t N, typename L> auto wholeOf(L &l)
  {
    using namespace WholeListDetail;
    using nakedL = std::remove_cvref_t<L>;
    auto fallback = [&] { return members<T>(l, std::make_index_sequence<N>()); };

    if constexpr (requires { l.wholeList(); })
      return l.wholeList();
    else if constexpr (requires { typename nakedL::IsWholeListSum; }) {
      // Sums and differences of two lists with equal member values.
      if constexpr (IsList<decltype(l.mR)> && IsList<decltype(l.mT)>) {
        auto a = wholeOf<T, NDim, N>(l.mR);
        auto b = wholeOf<T, NDim, N>(l.mT);
        using VA = MemberValue<decltype(a), NDim>;
        using VB = MemberValue<decltype(b), NDim>;
        constexpr int sign = nakedL::wholeListSign;
        if constexpr (std::is_same_v<VA, VB> || WidensTo<VB, VA>)
          return Sum<decltype(a), decltype(b)>{a, b, sign};
        else if constexpr (WidensTo<VA, VB>)
          // a + s b = s (b + s a), for s = +-1: accumulate in the wider type of b.
          return Scaled<double, Sum<decltype(b), decltype(a)>>{double(sign), {b, a, sign}};
        else
          return fallback();
      } else
        return fallback();
    } else if constexpr (requires { typename nakedL::IsWholeListProduct; }) {
      // A real scalar times a list.
      constexpr bool leftScalar = !IsList<decltype(l.mR)> && IsList<decltype(l.mT)>;
      constexpr bool rightScalar = IsList<decltype(l.mR)> && !IsList<decltype(l.mT)>;
      if constexpr (leftScalar || rightScalar) {
        auto &s = [&]() -> auto & {
          if constexpr (leftScalar)
            return l.mR;
          else
            return l.mT;
        }();
        auto &list = [&]() -> auto & {
          if constexpr (leftScalar)
            return l.mT;
          else
            return l.mR;
        }();
        using SV = ScalarValue<std::remove_cvref_t<decltype(s)>, NDim>;
        if constexpr (std::is_arithmetic_v<SV>) {
          auto a = wholeOf<T, NDim, N>(list);
          return Scaled<std::remove_cvref_t<decltype(s)>, decltype(a)>{s, a};
        } else
          return fallback();
      } else
        return fallback();
    } else if constexpr (requires { typename nakedL::IsWholeListNegation; }) {
      auto a = wholeOf<T, NDim, N>(l.mR);
      return Scaled<double, decltype(a)>{-1., a};
    } else
      return fallback();
  }
} // namespace TempLat

#endif
