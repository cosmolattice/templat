#ifndef TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_NUMBERCOLLECTION_H
#define TEMPLAT_LATTICE_ALGEBRA_CONSTANTS_NUMBERCOLLECTION_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#include <array>
#include <type_traits>
#include "TempLat/lattice/algebra/constants/number.h"
#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getcomponent.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{
  /** @brief A lightweight fixed-size collection of Number<T> values.
   *
   * Mirrors FieldCollection's interface for the evolver's collection-level
   * operations (Tag<I> accessors, operator+=, list-algebra integration),
   * but stores Number<T> instead of Field<T,NDIM>.
   *
   * Unit test: ctest -R test-numbercollection
   **/
  template <typename T, int N>
  struct NumberCollection {
    std::array<Number<T>, N> data{};

    // --- List-algebra integration ---
    using Getter = GetComponent;
    static constexpr size_t size = N;

    // --- Component access (makes IsTempLatGettable) ---
    template <int M>
    Number<T> &getComp(Tag<M>) { return data[M]; }

    template <int M>
    const Number<T> &getComp(Tag<M>) const { return data[M]; }

    // --- Tag<I> accessor (matches FieldCollection interface) ---
    template <int M>
      requires(M >= 0 && M < N)
    Number<T> &operator()(Tag<M>) { return data[M]; }

    template <int M>
      requires(M >= 0 && M < N)
    const Number<T> &operator()(Tag<M>) const { return data[M]; }

    // --- Element-wise operator+= from list expressions ---
    template <typename R>
    NumberCollection &operator+=(R &&r)
    {
      for_in_range<0, N>([&](auto i) { data[int(i)] += GetComponent::get(r, i); });
      return *this;
    }

    NumberCollection &operator+=(ZeroType) { return *this; }

    // --- Element-wise assignment from list expressions ---
    template <typename R>
      requires(!std::is_same_v<std::decay_t<R>, NumberCollection>)
    NumberCollection &operator=(R &&r)
    {
      for_in_range<0, N>(
          [&](auto i) { data[int(i)].value = DoEval::eval(GetComponent::get(r, i), size_t{0}); });
      return *this;
    }

    NumberCollection &operator=(const NumberCollection &other)
    {
      for (int i = 0; i < N; ++i)
        data[i].value = other.data[i].value;
      return *this;
    }
  };
} // namespace TempLat

#endif
