#ifndef TEMPLAT_PARALLEL_KOKKOS_INTERNAL_H
#define TEMPLAT_PARALLEL_KOKKOS_INTERNAL_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include "TempLat/parallel/devices/kokkos/kokkos.h"

namespace TempLat
{

  // Need to forward-declare this - we cannot include the layout header here, as it would create a circular
  // dependency.
  template <size_t NDim> struct LayoutStruct;

  namespace device_kokkos
  {
    // ------------------------------------------------
    // Getting View types with stars
    // ------------------------------------------------

    template <size_t NDim, typename T> struct GetKokkosNDStarType {
      using type = typename GetKokkosNDStarType<NDim - 1, T>::type *;
    };
    template <typename T> struct GetKokkosNDStarType<1, T> {
      using type = T *;
    };

    // ------------------------------------------------
    // View types
    // ------------------------------------------------

    using DefaultLayout = Kokkos::LayoutRight;
    using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
    using DefaultHostExecutionSpace = Kokkos::DefaultHostExecutionSpace;

    // ------------------------------------------------
    // Getting the policies for iteration
    // ------------------------------------------------

    template <size_t NDim> auto getLocalKokkosPolicy(const LayoutStruct<NDim> &layout)
    {
      Kokkos::Array<uint64_t, NDim> start_iteration;
      Kokkos::Array<uint64_t, NDim> stop_iteration;
      const auto localSizes = layout.getSizesInMemory();
      const size_t nGhosts = layout.getNGhosts();

      // What's going on here: on GPU, it is beneficial to reverse the memory access pattern, for coalesced access.
      // However, we do not want to impose this on the level of the memory layouts. In particular, this would
      // require additional transpositions when going to Fourier space, which is not what we want. So we do the
      // transposition within the thread dispatch, if we are on a GPU. Otherwise, for optimal cached memory access
      // on CPU, we do not reverse the access pattern.
      for (int d = 0; d < (int)NDim; ++d) {
        const int _d = device_kokkos::reverse_access_pattern ? (int)NDim - 1 - d : d;
        start_iteration[_d] = nGhosts;
        stop_iteration[_d] = start_iteration[_d] + localSizes[d];
      }

      if constexpr (NDim == 1) {
        return Kokkos::RangePolicy<DefaultExecutionSpace>(start_iteration[0], stop_iteration[0]);
      } else {
        return Kokkos::MDRangePolicy<DefaultExecutionSpace, Kokkos::Rank<NDim>>(start_iteration, stop_iteration);
      }
    }

    /**
     * @brief The CPU counterpart of getLocalKokkosPolicy: a policy over the OUTER NDim-1 dimensions
     * only. The contiguous last dimension is left to KokkosNDLambdaWrapperInnerLoop, which walks it
     * in a plain loop that the vectorizer can see. Never used on GPU, where the full-rank MDRange
     * with the reversed access pattern is what gives coalesced access.
     */
    template <size_t NDim> auto getLocalKokkosOuterPolicy(const LayoutStruct<NDim> &layout)
    {
      static_assert(NDim >= 2, "The outer policy needs a dimension left over for the inner loop.");
      const auto localSizes = layout.getSizesInMemory();
      const size_t nGhosts = layout.getNGhosts();

      if constexpr (NDim == 2) {
        return Kokkos::RangePolicy<DefaultExecutionSpace>(nGhosts, nGhosts + localSizes[0]);
      } else {
        Kokkos::Array<uint64_t, NDim - 1> start_iteration;
        Kokkos::Array<uint64_t, NDim - 1> stop_iteration;
        for (int d = 0; d < (int)NDim - 1; ++d) {
          start_iteration[d] = nGhosts;
          stop_iteration[d] = nGhosts + localSizes[d];
        }
        return Kokkos::MDRangePolicy<DefaultExecutionSpace, Kokkos::Rank<NDim - 1>>(start_iteration, stop_iteration);
      }
    }

    template <size_t NDim, typename I>
    auto getLocalKokkosPolicy(const device_kokkos::array<I, NDim> &starts, const device_kokkos::array<I, NDim> &stops)
    {
      Kokkos::Array<uint64_t, NDim> start_iteration;
      Kokkos::Array<uint64_t, NDim> stop_iteration;

      // What's going on here: on GPU, it is beneficial to reverse the memory access pattern, for coalesced access.
      // However, we do not want to impose this on the level of the memory layouts. In particular, this would
      // require additional transpositions when going to Fourier space, which is not what we want. So we do the
      // transposition within the thread dispatch, if we are on a GPU. Otherwise, for optimal cached memory access
      // on CPU, we do not reverse the access pattern.
      for (int d = 0; d < (int)NDim; ++d) {
        const int _d = device_kokkos::reverse_access_pattern ? (int)NDim - 1 - d : d;
        start_iteration[_d] = starts[d];
        stop_iteration[_d] = start_iteration[_d] + stops[d];
      }

      if constexpr (NDim == 1) {
        return Kokkos::RangePolicy<DefaultExecutionSpace>(start_iteration[0], stop_iteration[0]);
      } else {
        return Kokkos::MDRangePolicy<DefaultExecutionSpace, Kokkos::Rank<NDim>>(start_iteration, stop_iteration);
      }
    }

    // ------------------------------------------------
    // Helpers for tuples and arrays
    // ------------------------------------------------

    // If you wonder about how performant this is, take a look at https://godbolt.org/z/d33cMaEG1
    /**
     * @brief Returns the last i elements of the given tuple t
     *
     * @tparam i size of the tuple to be returned
     * @param t tuple to be split
     * @return auto a tied tuple of the last i elements
     */
    template <size_t i, typename Head, typename... Tail>
      requires(i <= sizeof...(Tail))
    DEVICE_INLINE_FUNCTION constexpr auto tuple_last(const device_kokkos::tuple<Head, Tail...> &t)
    {
      static_assert(i <= sizeof...(Tail), "Cannot take a longer tail than the tuple.");
      if constexpr (sizeof...(Tail) + 1 == i)
        return apply([](auto &head, auto &...tail) { return device_kokkos::tie(head, tail...); }, t);
      else if constexpr (sizeof...(Tail) == i)
        return apply([](auto & /*head*/, auto &...tail) { return device_kokkos::tie(tail...); }, t);
      else
        return apply(
            [](auto & /*head*/, auto &...tail) { return device_kokkos::tuple_last<i>(device_kokkos::tie(tail...)); },
            t);
    }

    /**
     * @brief Returns the first i elements of the given tuple t
     *
     * @tparam i size of the tuple to be returned
     * @param t tuple to be split
     * @return auto a tied tuple of the first i elements
     */
    template <int i, typename Head, typename... Tail>
    DEVICE_INLINE_FUNCTION constexpr auto tuple_first(const device_kokkos::tuple<Head, Tail...> &t)
    {
      static_assert(i <= sizeof...(Tail), "Cannot take a longer sequence than the tuple.");
      static_assert(i >= 0, "Cannot take a longer sequence than the tuple.");
      if constexpr (i == 0)
        return device_kokkos::tuple();
      else if constexpr (i == 1)
        return apply([](auto &head, auto &.../*tail*/) { return device_kokkos::tie(head); }, t);
      else
        return apply(
            [](auto &head, auto &...tail) {
              return device_kokkos::tuple_cat(device_kokkos::tie(head),
                                              device_kokkos::tuple_first<i - 1>(device_kokkos::tie(tail...)));
            },
            t);
    }

    template <typename... Args, std::size_t... Is>
    DEVICE_INLINE_FUNCTION auto reverse_tuple(const device_kokkos::tuple<Args...> &tuple, std::index_sequence<Is...>)
    {
      return device_kokkos::tie(device_kokkos::get<sizeof...(Args) - 1 - Is>(tuple)...);
    }

    template <typename... Args> DEVICE_INLINE_FUNCTION auto reverse_tuple(const device_kokkos::tuple<Args...> &tuple)
    {
      return reverse_tuple(tuple, std::make_index_sequence<sizeof...(Args)>());
    }

    template <typename Arg, size_t N, std::size_t... Is>
    DEVICE_INLINE_FUNCTION auto reverse_array(const device_kokkos::array<Arg, N> &arr, std::index_sequence<Is...>)
    {
      return device_kokkos::array<Arg, N>{{get<N - 1 - Is>(arr)...}};
    }

    template <typename Arg, size_t N> DEVICE_INLINE_FUNCTION auto reverse_array(const device_kokkos::array<Arg, N> &arr)
    {
      return reverse_array(arr, std::make_index_sequence<N>());
    }

    // ------------------------------------------------
    // Wrapping multidimensional device lambdas for use with Kokkos
    // ------------------------------------------------

    /**
     * @brief This is a functor which wraps a lambda.
     * Basically, this is necessary when one wants to call a variadic lambda on an NVIDIA GPU.
     * CUDA seems to be unable to expand the variadic arguments - in contrast, a direct approach does indeed work for
     * openMP or serial compilation.
     * To get around this limitation, the KokkosNDLambdaWrapper packs the indices into an array.
     * If you wonder, whether there's a difference when using tie and tuples: https://godbolt.org/z/M3bG39rsM
     * No. Therefore, we spare the ourselves the hassle and simply use an array.
     *
     * @tparam NDim Number of arguments taken
     * @tparam FUN The lambda to which we forward the indices
     */
    template <size_t NDim, typename FUN> struct KokkosNDLambdaWrapper {
      KokkosNDLambdaWrapper(const FUN &_fun) : fun(_fun) {};

      template <typename... Args>
        requires(sizeof...(Args) == NDim)
      DEVICE_INLINE_FUNCTION void operator()(const Args &...args) const
      {
        // What's going on here: on GPU, it is beneficial to reverse the memory access pattern, for coalesced access.
        // However, we do not want to impose this on the level of the memory layouts. In particular, this would
        // require additional transpositions when going to Fourier space, which is not what we want. So we do the
        // transposition within the thread dispatch, if we are on a GPU. Otherwise, for optimal cached memory access
        // on CPU, we do not reverse the access pattern.
        if constexpr (device_kokkos::reverse_access_pattern)
          fun(device_kokkos::reverse_array(device_kokkos::IdxArray<NDim>{{static_cast<Idx>(args)...}}));
        else
          fun({{static_cast<Idx>(args)...}});
      }

      FUN fun;
    };

    /**
     * @brief The CPU dispatch: a functor which wraps a lambda and walks the contiguous last
     * dimension itself.
     *
     * Kokkos' MDRange tiles all NDim dimensions, and its innermost tile loop is one the vectorizer
     * gives up on -- it has runtime bounds and no independence guarantee, so GCC emits scalar code
     * for every lattice kernel. Handing the last (unit-stride) dimension back to a plain loop that
     * carries TEMPLAT_ASSUME_INDEPENDENT is what lets it vectorize: measured on the SU(2) kick, 125
     * -> 46 ns/site, and 8 -> 8064 packed FP instructions. See benchmarks/PERFORMANCE.md.
     *
     * The functor is copied once, here at construction; the inner loop must not copy it again. The
     * expression trees it holds are NOT cheap to copy (their leaves carry a shared_ptr to the memory
     * manager, so a copy means atomic refcount traffic), and copying one per row costs more than
     * vectorising it saves.
     *
     * @tparam NDim Dimension of the lattice; the wrapped functor still takes all NDim indices.
     * @tparam FUN The lambda to which we forward the indices
     */
    template <size_t NDim, typename FUN> struct KokkosNDLambdaWrapperInnerLoop {
      KokkosNDLambdaWrapperInnerLoop(const FUN &_fun, Idx _innerStart, Idx _innerStop)
          : fun(_fun), innerStart(_innerStart), innerStop(_innerStop) {};

      template <typename... Args>
        requires(sizeof...(Args) == NDim - 1)
      DEVICE_INLINE_FUNCTION void operator()(const Args &...args) const
      {
        const Idx start = innerStart, stop = innerStop;
        TEMPLAT_ASSUME_INDEPENDENT
        for (Idx i = start; i < stop; ++i)
          fun(device_kokkos::IdxArray<NDim>{{static_cast<Idx>(args)..., i}});
      }

      FUN fun;
      Idx innerStart, innerStop;
    };

    /**
     * @brief This is a functor which wraps a lambda for reduction.
     * Basically, this is necessary when one wants to call a variadic lambda on an NVIDIA GPU.
     * CUDA seems to be unable to expand the variadic arguments - in contrast, a direct approach does indeed work for
     * openMP or serial compilation.
     * To get around this limitation, the KokkosNDLambdaWrapperReduction packs the indices into an array.
     * If you wonder, whether there's a difference when using tie and tuples: https://godbolt.org/z/M3bG39rsM
     * No. Therefore, we spare the ourselves the hassle and simply use an array.
     *
     * @tparam NDim Number of arguments taken
     * @tparam FUN The lambda to which we forward the indices
     */
    template <size_t NDim, typename FUN> struct KokkosNDLambdaWrapperReduction {
      KokkosNDLambdaWrapperReduction(const FUN &_fun) : fun(_fun) {};

      template <typename... Args>
        requires(sizeof...(Args) == NDim + 1)
      DEVICE_INLINE_FUNCTION void operator()(Args &&...args) const
      {
        // What's going on here: on GPU, it is beneficial to reverse the memory access pattern, for coalesced access.
        // However, we do not want to impose this on the level of the memory layouts. In particular, this would
        // require additional transpositions when going to Fourier space, which is not what we want. So we do the
        // transposition within the thread dispatch, if we are on a GPU. Otherwise, for optimal cached memory access
        // on CPU, we do not reverse the access pattern.
        auto tuple = device_kokkos::tie(args...);
        if constexpr (device_kokkos::reverse_access_pattern)
          fun(makeArray(device_kokkos::reverse_tuple(device_kokkos::tuple_first<NDim>(tuple))),
              device_kokkos::get<NDim>(tuple)); // the last argument is the reduction result
        else
          fun(makeArray(device_kokkos::tuple_first<NDim>(tuple)),
              device_kokkos::get<NDim>(tuple)); // the last argument is the reduction result
      }

      FUN fun;

      template <typename... Args>
        requires(sizeof...(Args) == NDim)
      DEVICE_INLINE_FUNCTION auto makeArray(device_kokkos::tuple<Args...> &&tuple) const
      {
        return device_kokkos::apply(
            [](const auto &...args) { return device_kokkos::IdxArray<NDim>{{static_cast<Idx>(args)...}}; },
            std::move(tuple));
      }
    };
  } // namespace device_kokkos
} // namespace TempLat

#endif // KOKKOS_H
