#ifndef TEMPLAT_PARALLEL_KOKKOS_MEMORY_H
#define TEMPLAT_PARALLEL_KOKKOS_MEMORY_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"

#include "TempLat/parallel/devices/kokkos/kokkos.h"
#include "TempLat/parallel/devices/kokkos/kokkos_internal.h"

#include <Kokkos_Core.hpp>
#include <sstream>

namespace TempLat::device_kokkos::memory
{
  template <typename T, size_t NDim, typename Exec = DefaultExecutionSpace, typename Layout = DefaultLayout>
  using NDView = Kokkos::View<typename GetKokkosNDStarType<NDim, T>::type, // Get the star syntax for
                                                                           // dimensionality recursively with
                              Layout, // LayoutRight is most compatible for now, may change in future
                              Exec    // Choice between GPU and CPU
                              >;
  template <typename T, size_t NDim, typename Exec = DefaultExecutionSpace, typename Layout = DefaultLayout>
  using NDViewUnmanaged =
      Kokkos::View<typename GetKokkosNDStarType<NDim, T>::type, // Get the star syntax for dimensionality
                                                                // recursively
                   Layout, // LayoutRight is most compatible for now, may change in future
                   Exec,   // Choice between GPU and CPU
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> // No allocation: Attach to existing memory
                   >;

  template <typename T, size_t NDim>
  using NDViewUnmanagedHost =
      Kokkos::View<typename GetKokkosNDStarType<NDim, T>::type, // Get the star syntax for dimensionality
                                                                // recursively
                   DefaultLayout,             // LayoutRight is most compatible for now, may change in future
                   DefaultHostExecutionSpace, // Choice between GPU and CPU
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> // No allocation: Attach to existing memory
                   >;

  using Kokkos::subview;

  template <typename A> auto createMirrorView(const A &a) { return Kokkos::create_mirror_view(a); }

  template <typename OBJ, size_t NDim, typename T>
  void setAtOnePoint(OBJ &&obj, device_kokkos::IdxArray<NDim> pos, T val)
  {
    Kokkos::parallel_for(
        "Set a point", Kokkos::RangePolicy(0, 1), DEVICE_LAMBDA(const unsigned int) {
          device_kokkos::apply([&](const auto... idx) { obj.getSet(idx...) = val; }, pos);
        });
  }

  template <typename View, typename T> void fill(View &view, const T &value) { Kokkos::deep_copy(view, value); }

  template <typename OBJ, size_t NDim, typename I = ptrdiff_t>
  GetGetReturnType<OBJ>::type getAtOnePoint(OBJ &&obj, const device_kokkos::array<I, NDim> &pos)
  {
    using T = GetGetReturnType<OBJ>::type;
    T ret;
    Kokkos::parallel_reduce(
        "Get a point", Kokkos::RangePolicy(0, 1),
        DEVICE_LAMBDA(const unsigned int, T &update) {
          device_kokkos::apply([&](const auto... idx) { update = DoEval::eval(obj, idx...); }, pos);
        },
        ret);
    return ret;
  }

  template <typename View1, typename View2>
    requires(Kokkos::is_view<View1>::value && Kokkos::is_view<View2>::value)
  void copyDeviceToDevice(const View1 &src, View2 &dest)
  {
    static_assert(View1::rank == View2::rank, "Source and destination views must have the same rank.");
    static_assert(std::is_same_v<typename View1::value_type, typename View2::value_type>,
                  "Source and destination views must have the same value type.");
    static_assert(std::is_same_v<typename View1::execution_space, typename View2::execution_space>,
                  "Source and destination views must have the same execution space.");
    // static_assert(std::is_same_v<typename View1::array_layout, typename View2::array_layout>,
    //               "Source and destination views must have the same layout.");

    constexpr size_t dim = View1::rank;
    for (size_t i = 0; i < dim; ++i)
      if (src.extent(i) != dest.extent(i)) {
        std::stringstream ss;
        ss << "Source and destination views must have the same extents. Mismatch at dimension " << i << ": "
           << "src extent = " << src.extent(i) << ", dest extent = " << dest.extent(i);
        throw std::runtime_error(ss.str());
      }

    bool contiguous = src.span_is_contiguous() && dest.span_is_contiguous();
    if (contiguous) {
      // Stream-ordered copy instead of a fenced one.
      //
      // The no-instance Kokkos::deep_copy overload wraps the copy in two DEVICE-WIDE
      // Kokkos::fence() calls (Kokkos 5.1.1, Kokkos_CopyViews.hpp:1316-1325 for the
      // byte-wise path, 1327-1331 for the view_copy fallback). Ghost exchange calls
      // this ~400x per MC sweep -- two of the six face copies per component per
      // parity pass land here, the ones whose slab happens to be contiguous -- which
      // measured as ~800 cudaDeviceSynchronize/sweep and a large part of the ~14 ms
      // volume-independent floor in the sweep cost.
      //
      // The instance overload does not fence on either path for our case: the
      // byte-wise path (2517-2521) just enqueues DeepCopy on the instance, and the
      // fallback takes Impl::view_copy(exec_space, ...) because Cuda can access
      // CudaSpace on both sides.
      //
      // Ordering w.r.t. the sweep kernels holds because TempLat only ever uses the
      // default execution space instance, i.e. one CUDA stream -- the non-contiguous
      // branch below is a bare parallel_for and has always depended on exactly this.
      // Callers needing host- or MPI-visible completion must fence themselves, and
      // the ghost pack/unpack already does (ghostupdater.h:321,337,446,465); host
      // reads go through copyDeviceToHost, whose deep_copy fences on its own.
      //
      // Revisit if a second execution space instance / stream is ever introduced.
      Kokkos::deep_copy(typename View1::execution_space(), dest, src);
    } else {
      // If not, we need to do a manual copy
      device::array<ptrdiff_t, dim> localSizes;
      for (size_t i = 0; i < dim; ++i)
        localSizes[i] = src.extent(i);

      auto functor = DEVICE_LAMBDA(const device_kokkos::IdxArray<dim> &idx)
      {
        device_kokkos::apply([&](const auto... i) { dest(i...) = src(i...); }, idx);
      };

      Kokkos::parallel_for("Copy non-contiguous", getLocalKokkosPolicy({}, localSizes),
                           KokkosNDLambdaWrapper<dim, decltype(functor)>(functor));
    }
  }

  template <typename View>
    requires(Kokkos::is_view<View>::value)
  void copyDeviceToHost(View &src, typename View::value_type *dest)
  {
    constexpr size_t dim = View::rank();
    using T = typename View::value_type;
    using Exec = typename View::execution_space;
    using Layout = typename View::array_layout;

    const bool contiguous = src.span_is_contiguous();

    // If the source view is contiguous, we can use a simple copy
    if (contiguous) {
      auto destView = NDViewUnmanaged<T, dim, Exec, Layout>(dest, src.layout());
      Kokkos::deep_copy(destView, src);
    } else {
      // If not, we first need a temporary contiguous copy on device
      device_kokkos::array<ptrdiff_t, dim> localSizes;
      for (size_t i = 0; i < dim; ++i)
        localSizes[i] = src.extent(i);
      auto device_temp = device_kokkos::apply(
          [&](const auto... sizes) { return NDView<T, dim, Exec, Layout>("temp", sizes...); }, localSizes);
      copyDeviceToDevice(src, device_temp);

      // now copy the contiguous temp to host
      auto destView = NDViewUnmanaged<T, dim, Exec, Layout>(dest, device_temp.layout());
      Kokkos::deep_copy(destView, device_temp);
    }
  }

  template <typename View>
    requires(Kokkos::is_view<View>::value)
  void copyHostToDevice(const typename View::value_type *src, View &dest)
  {
    constexpr size_t dim = View::rank();
    using T = typename View::value_type;
    using Exec = typename View::execution_space;
    using Layout = typename View::array_layout;
    const bool contiguous = dest.span_is_contiguous();

    // If the destination view is contiguous, we can use a simple copy
    if (contiguous) {
      auto srcView = NDViewUnmanaged<T, dim, Exec, Layout>(const_cast<T *>(src), dest.layout());
      Kokkos::deep_copy(dest, srcView);
    } else {
      // If not, we first need a temporary contiguous copy on device
      device_kokkos::array<ptrdiff_t, dim> localSizes;
      for (size_t i = 0; i < dim; ++i)
        localSizes[i] = dest.extent(i);
      auto device_temp = device_kokkos::apply(
          [&](const auto... sizes) { return NDView<T, dim, Exec, Layout>("temp", sizes...); }, localSizes);
      copyHostToDevice(src, device_temp);

      // now copy the contiguous temp to the original view
      copyDeviceToDevice(device_temp, dest);
    }
  }
} // namespace TempLat::device_kokkos::memory

#endif
