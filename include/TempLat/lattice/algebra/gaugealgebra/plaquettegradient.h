#ifndef TEMPLAT_LATTICE_ALGEBRA_GAUGEALGEBRA_PLAQUETTEGRADIENT_H
#define TEMPLAT_LATTICE_ALGEBRA_GAUGEALGEBRA_PLAQUETTEGRADIENT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <string>
#include <type_traits>
#include <utility>

#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquette.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquetteback.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getcomponent.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/parallel/device.h"
#include "TempLat/util/constexpr_for.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/staticif.h"

namespace TempLat
{
  namespace PlaqGradDetail
  {
    /** @brief The p-th pair (i, j), 1 <= i < j <= NDim, in lexicographic order. */
    template <size_t NDim> constexpr std::pair<int, int> upperPair(size_t p)
    {
      for (int i = 1; i <= (int)NDim; ++i)
        for (int j = i + 1; j <= (int)NDim; ++j)
          if (p-- == 0) return {i, j};
      return {0, 0};
    }
    /** @brief The q-th pair (i, j), i != j, in lexicographic order. */
    template <size_t NDim> constexpr std::pair<int, int> orderedPair(size_t q)
    {
      for (int i = 1; i <= (int)NDim; ++i)
        for (int j = 1; j <= (int)NDim; ++j)
          if (i != j && q-- == 0) return {i, j};
      return {0, 0};
    }

    /** @brief plaqGrad evaluated for all directions at once, see PlaqGrad::wholeList. */
    template <size_t NDim, typename Forward, typename Backward> struct Whole {
      static constexpr size_t size = NDim;
      static constexpr bool shares = true;
      Forward forward;   // plaq(Us, i, j) for i < j, in upperPair order
      Backward backward; // plaqBack(Us, i, j) for i != j, in orderedPair order

      template <typename... IDX> DEVICE_INLINE_FUNCTION auto evalWhole(const IDX &...idx) const
      {
        using V = std::remove_cvref_t<decltype(DoEval::eval(device::get<0>(forward), idx...))>;
        device::array<V, NDim> out{};
        // P_ji = P_ij^dagger: the components after the first (the identity) change sign, for SU(2) and U(1).
        constexpr_for<0, NDim *(NDim - 1) / 2>([&](auto p) {
          constexpr auto ij = upperPair<NDim>(decltype(p)::value);
          const V v = DoEval::eval(device::get<decltype(p)::value>(forward), idx...);
          for (size_t c = 0; c < v.size(); ++c) {
            out[ij.first - 1][c] += v[c];
            out[ij.second - 1][c] += c == 0 ? v[c] : -v[c];
          }
        });
        constexpr_for<0, NDim *(NDim - 1)>([&](auto q) {
          constexpr auto ij = orderedPair<NDim>(decltype(q)::value);
          const V v = DoEval::eval(device::get<decltype(q)::value>(backward), idx...);
          for (size_t c = 0; c < v.size(); ++c)
            out[ij.first - 1][c] -= v[c];
        });
        return out;
      }

      void preGet()
      {
        device::apply([](auto &...f) { (PreGet::apply(f), ...); }, forward);
        device::apply([](auto &...b) { (PreGet::apply(b), ...); }, backward);
      }
      void postGet()
      {
        device::apply([](auto &...f) { (PostGet::apply(f), ...); }, forward);
        device::apply([](auto &...b) { (PostGet::apply(b), ...); }, backward);
      }
    };
  } // namespace PlaqGradDetail

  /** @brief The list over directions i = 1 ... NDim whose member i is
   *  sum_{j != i} [plaq(Us, i, j) - plaqBack(Us, i, j)]: the sum of the plaquettes containing the link U_i(x), each
   *  starting with it (plaqBack runs the other way). Up to normalization the gradient of the plaquette action with
   *  respect to U_i, i.e. the lattice D_j F_ji; projecting onto the algebra is left to the assignment.
   *
   *  Member by member it is the expression above. A collection assignment evaluates it for all directions at once
   *  (wholeList), computing every forward plaquette once, since P_ji = P_ij^dagger: in three dimensions 9 instead
   *  of 12 four-link loops per site.
   *
   *  Unit test: ctest -R test-plaquettegradient
   **/
  template <typename R> class PlaqGrad
  {
  public:
    static constexpr size_t NDim = R::NDim;
    static constexpr size_t size = NDim;
    using Getter = GetComponent;

    explicit PlaqGrad(const R &Us) : mUs(Us) {}

    template <int M> auto getComp(Tag<M>) const
    {
      constexpr Tag<M + 1> i;
      return Total(j, 1, (int)NDim, IfElse(i != j, plaq(mUs, i, j) - plaqBack(mUs, i, j), ZeroType()));
    }

    /** @brief The evaluation of all directions at once, for collection assignments. */
    auto wholeList() const
    {
      using namespace PlaqGradDetail;
      auto forward = [&]<size_t... p>(std::index_sequence<p...>) {
        return device::make_tuple(plaq(mUs, Tag<upperPair<NDim>(p).first>(), Tag<upperPair<NDim>(p).second>())...);
      }(std::make_index_sequence<NDim *(NDim - 1) / 2>());
      auto backward = [&]<size_t... q>(std::index_sequence<q...>) {
        return device::make_tuple(
            plaqBack(mUs, Tag<orderedPair<NDim>(q).first>(), Tag<orderedPair<NDim>(q).second>())...);
      }(std::make_index_sequence<NDim *(NDim - 1)>());
      return Whole<NDim, decltype(forward), decltype(backward)>{forward, backward};
    }

    std::string toString() const { return "plaqGrad"; }

  private:
    R mUs;
  };

  /**
   * @vocab-summary The list over directions $i$ of $\sum_{j\neq i}[P_{ij}(x) - P^{\rm back}_{ij}(x)]$, the
   * plaquettes containing $U_i(x)$: up to normalization the gradient of the plaquette action, i.e. the lattice
   * $D_j F_{ji}$. Assigned to a collection, each forward plaquette is computed once.
   * @vocab-signature plaqGrad(Us)
   **/
  template <typename R> PlaqGrad<R> plaqGrad(const R &Us) { return PlaqGrad<R>(Us); }
} // namespace TempLat

#endif
