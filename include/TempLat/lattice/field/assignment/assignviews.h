#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNMENT_ASSIGNVIEWS_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNMENT_ASSIGNVIEWS_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <algorithm>
#include <string>
#include <type_traits>

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/ghostreach.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstruct.h"
#include "TempLat/lattice/memory/sweepbuffer.h"
#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/util/constexpr_for.h"

namespace TempLat
{
  /** @brief Number of planes evaluated per sweep step: about 8 sites per resident device thread, at least one
   *  plane and at most all of them. A positive requested value overrides the choice.
   *
   *  The factor 8 comes from bench-sweep_assign (256^3, RTX 4070): heavy stencils keep fewer threads resident
   *  than the device maximum, and the sweep time levels off only from about 8 sites per thread. On a CPU this
   *  gives one plane. */
  inline device::Idx sweepBandWidth(device::Idx planeArea, device::Idx nPlanes, device::Idx requested)
  {
    constexpr device::Idx waves = 8;
    const device::Idx sites = waves * device::iteration::concurrency();
    const device::Idx width = requested > 0 ? requested : (sites + planeArea - 1) / planeArea;
    return std::clamp<device::Idx>(width, 1, nPlanes);
  }

  /** @brief A right-hand side evaluating many fields per site: a whole collection, see ComponentList. On a CPU, its
   *  kernel inlines everything it calls; GCC otherwise calls the member expressions per site and recomputes the
   *  index arithmetic of every field access. */
  template <typename R>
  concept FusesComponents = requires { requires std::remove_cvref_t<R>::fusesComponents; };

  /** @brief Whether the kernel assigning g to NComp views stores member by member: a collection of composite fields,
   *  each member filling NComp / size views. That keeps a single member's intermediates live at a time, which on a
   *  CPU removes most register spills of e.g. a fused gauge kick. A collection of plain fields stores all members
   *  at once. */
  template <typename R, size_t NComp> constexpr bool storesPerMember()
  {
    if constexpr (FusesComponents<R>)
      return NComp > std::remove_cvref_t<R>::size;
    else
      return false;
  }

#if defined(__GNUC__) && !defined(__CUDACC__) && !defined(__HIPCC__)
#define TEMPLAT_INLINE_CALLEES [[gnu::flatten]]
#else
#define TEMPLAT_INLINE_CALLEES
#endif

  namespace AssignViewsDetail
  {
    template <typename View, size_t NDim>
    DEVICE_INLINE_FUNCTION decltype(auto) at(const View &v, const device::IdxArray<NDim> &i)
    {
      return device::apply([&](auto... k) -> decltype(auto) { return v(k...); }, i);
    }

    /** @brief The assignment kernel: stores component first + c of g's value at a site into views[c]; a scalar
     *  value goes to views[0]. For a collection of composite fields (storesPerMember), member m fills
     *  views[m * K + k] with component first + k of its value. */
    template <size_t first, typename View, size_t NComp, typename R> struct AssignKernel {
      static constexpr size_t NDim = View::rank;
      /** @brief Asks the iteration to inline this kernel into its loop, see KokkosNDLambdaWrapper. A collection of
       *  composite fields needs it; for plain fields GCC inlines it anyway, and forcing costs ~4 % (fused N = 8 kick). */
      static constexpr bool inlineIntoLoop = storesPerMember<R, NComp>();
      device::array<View, NComp> views;
      R g;

      DEVICE_INLINE_FUNCTION void operator()(const device::IdxArray<NDim> &idx) const
        requires(!FusesComponents<R>)
      {
        device::apply([&](const auto &...i) { store(i...); }, idx);
      }
      TEMPLAT_INLINE_CALLEES DEVICE_FORCEINLINE_FUNCTION void operator()(const device::IdxArray<NDim> &idx) const
        requires FusesComponents<R>
      {
        device::apply([&](const auto &...i) { store(i...); }, idx);
      }

      template <typename... IDX> DEVICE_FORCEINLINE_FUNCTION void store(const IDX &...i) const
      {
        if constexpr (storesPerMember<R, NComp>()) {
          constexpr size_t K = NComp / R::size;
          constexpr_for<0, R::size>([&](auto m) {
            const auto value = g.template evalMember<decltype(m)::value>(i...);
            constexpr_for<0, K>([&](auto k) { views[decltype(m)::value * K + k](i...) = value[first + k]; });
          });
        } else {
          const auto value = DoEval::eval(g, i...);
          if constexpr (std::is_convertible_v<decltype(value), typename View::value_type>) {
            static_assert(NComp == 1 && first == 0, "A scalar value fills a single view.");
            views[0](i...) = value;
          } else
            constexpr_for<0, NComp>([&](auto c) { views[c](i...) = value[first + c]; });
        }
      }
    };

    /** @brief Copies sites from the ring of a sweep back into the views: plane p (counted from the first local
     *  plane) sits in ring plane p mod nSlots. */
    template <typename View, size_t NComp> struct WriteBackKernel {
      static constexpr size_t NDim = View::rank;
      device::array<View, NComp> views, ring;
      device::Idx nSlots, nGhosts;

      DEVICE_INLINE_FUNCTION void operator()(const device::IdxArray<NDim> &idx) const
      {
        device::IdxArray<NDim> slot = idx;
        slot[0] = (idx[0] - nGhosts) % nSlots;
        for (size_t c = 0; c < NComp; ++c)
          at(views[c], idx) = at(ring[c], slot);
      }
    };

    /** @brief A view of the same type and shape as v, on the memory at data. */
    template <typename View> View sameShape(const View &v, typename View::value_type *data)
    {
      device::IdxArray<View::rank> extents;
      for (size_t d = 0; d < View::rank; ++d)
        extents[d] = v.extent(d);
      return device::apply([&](auto... n) { return View(data, n...); }, extents);
    }

    /** @brief The bands an assignment is evaluated in, see assignViews: either one band covering the local block,
     *  written directly into the views, or a sweep through a ring buffer. Independent of the right-hand side. */
    template <typename View, size_t NComp> class Bands
    {
    public:
      static constexpr size_t NDim = View::rank;
      using Views = device::array<View, NComp>;

      /** @brief A single band. */
      Bands(const Views &views, const LayoutStruct<NDim> &layout)
          : mViews(views), mLayout(layout), mNPlanes(layout.getSizesInMemory()[0]), mBand(mNPlanes), mNBands(1)
      {
      }

      /** @brief A sweep for a right-hand side reading the views down to offset -reachDown along axis 0. */
      Bands(const Views &views, const LayoutStruct<NDim> &layout, device::Idx reachDown, SweepBuffer &buffer)
          : mViews(views), mLayout(layout), mNPlanes(layout.getSizesInMemory()[0]), mSweep(true)
      {
        using device::Idx;
        const auto sizes = layout.getSizesInMemory();
        Idx planeArea = 1;
        for (size_t d = 1; d < NDim; ++d)
          planeArea *= sizes[d];
        mBand = sweepBandWidth(planeArea, mNPlanes, buffer.getBandWidth());
        mNBands = (mNPlanes + mBand - 1) / mBand;
        // Band b occupies ring block b mod nBlocks, i.e. ring planes (b mod nBlocks) * band + [0, band), so
        // plane p sits in ring plane p mod nSlots. A block is written back just before band b + nBlocks reuses
        // it; that band reads none of its planes.
        mNBlocks = std::min<Idx>(1 + (reachDown + mBand - 1) / mBand, mNBands);
        const Idx nSlots = mNBlocks * mBand;

        mPlaneStride = views[0].stride(0);
        auto *ringData = buffer.get<typename View::value_type>(NComp * nSlots * mPlaneStride);
        Views ring;
        for (size_t c = 0; c < NComp; ++c) {
          device::IdxArray<NDim> extents;
          extents[0] = nSlots;
          for (size_t d = 1; d < NDim; ++d)
            extents[d] = views[c].extent(d);
          ring[c] = device::apply([&](auto... n) { return View(ringData + c * nSlots * mPlaneStride, n...); }, extents);
        }
        mWriteBack = {views, ring, nSlots, layout.getNGhosts()};
      }

      device::Idx count() const { return mNBands; }
      device::Idx firstPlane(device::Idx b) const { return b * mBand; }
      device::Idx planes(device::Idx b) const { return std::min(mBand, mNPlanes - firstPlane(b)); }

      /** @brief Where band b is written: the views, or views of the same shape shifted into the ring. Before
       *  band b reuses ring planes, writes their content back. */
      Views prepare(device::Idx b)
      {
        if (!mSweep) return mViews;
        if (b >= mNBlocks) {
          device::iteration::foreach ("SweepWriteBack", mLayout, mWritten, mBand, mWriteBack);
          mWritten += mBand;
        }
        // Memory plane nGhosts + firstPlane(b) maps to ring plane (b mod nBlocks) * band.
        const device::Idx shift = (b % mNBlocks) * mBand - (mLayout.getNGhosts() + firstPlane(b));
        Views target;
        for (size_t c = 0; c < NComp; ++c)
          target[c] = sameShape(mViews[c], mWriteBack.ring[c].data() + shift * mPlaneStride);
        return target;
      }

      /** @brief Writes back the planes still in the ring. */
      void finish()
      {
        if (mSweep) device::iteration::foreach ("SweepWriteBack", mLayout, mWritten, mNPlanes - mWritten, mWriteBack);
      }

    private:
      Views mViews;
      const LayoutStruct<NDim> &mLayout;
      device::Idx mNPlanes, mBand, mNBands, mNBlocks = 1, mPlaneStride = 0, mWritten = 0;
      bool mSweep = false;
      WriteBackKernel<View, NComp> mWriteBack{};
    };
  } // namespace AssignViewsDetail

  /** @brief Assigns g to the component views of a field, as one kernel: views[c] receives component first + c of
   *  g's value, a scalar value goes to views[0]. targetReach is the range of offsets at which g reads the views,
   *  see prepareAssignment.
   *
   *  If g reads the views only in place, the kernel runs once over the local block. Otherwise, or if a kernel
   *  storing member by member has a member reading another member's views (readsOtherMember, see
   *  prepareCollectionAssignment), it sweeps along
   *  axis 0, the slowest in memory, in bands of B adjacent planes. The kernel evaluates each band into a ring
   *  buffer, never into the views, so every read sees old values: it is handed views of the same type and shape
   *  whose data pointer is shifted to put the band into the ring. Before a band reuses ring planes, a separate
   *  kernel writes their content back; that happens only once they lie below everything the band reads.
   *  Reads above the band hit untouched planes, reads past the local block land on ghosts, which are never
   *  written.
   *
   *  The plain and the sweeping assignment share the kernel and its launch site; everything else does not depend
   *  on g. Extra memory: NComp * B * (1 + ceil(r / B)) planes, r being the downward reach along axis 0.
   *  Launches: two per band.
   *
   *  Unit test: ctest -R test-sweepassignment
   **/
  template <size_t first = 0, typename View, size_t NComp, typename R, size_t NDim>
  void assignViews(const std::string &name, const device::array<View, NComp> &views, R &&g,
                   const LayoutStruct<NDim> &layout, const GhostReach &targetReach, SweepBuffer &buffer,
                   bool readsOtherMember = false)
  {
    using Kernel = AssignViewsDetail::AssignKernel<first, View, NComp, std::remove_cvref_t<R>>;

    PreGet::apply(g);

    // One launch site, so that the plain and the sweeping assignment share all code depending on g.
    using Bands = AssignViewsDetail::Bands<View, NComp>;
    const bool sweep = !targetReach.isZero() || (storesPerMember<R, NComp>() && readsOtherMember);
    Bands bands = sweep ? Bands(views, layout, -targetReach.lo[0], buffer) : Bands(views, layout);
    // Row tiles help the many streams of a collection of plain fields; stencils of composites prefer Kokkos' tiles.
    constexpr bool rowTiles = FusesComponents<R> && !storesPerMember<R, NComp>();
    for (device::Idx b = 0; b < bands.count(); ++b)
      device::iteration::foreach (name, layout, bands.firstPlane(b), bands.planes(b), Kernel{bands.prepare(b), g},
                                  rowTiles);
    bands.finish();

    PostGet::apply(g);
  }
} // namespace TempLat

#endif
