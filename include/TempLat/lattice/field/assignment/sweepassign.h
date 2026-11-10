#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNMENT_SWEEPASSIGN_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNMENT_SWEEPASSIGN_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <algorithm>

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstruct.h"
#include "TempLat/lattice/memory/sweepbuffer.h"
#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"

namespace TempLat
{
  namespace SweepAssignDetail
  {
    template <typename View, size_t NDim>
    DEVICE_INLINE_FUNCTION decltype(auto) at(const View &v, const device::IdxArray<NDim> &i)
    {
      return device::apply([&](auto... k) -> decltype(auto) { return v(k...); }, i);
    }
  } // namespace SweepAssignDetail

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

  /** @brief Assigns g to view in place, where g may read view at neighbouring sites.
   *
   *  Sweeps along axis 0, the slowest in memory. Each step evaluates a band of adjacent planes into a ring
   *  buffer, never into view, so every read of view sees old values. A ring slot is reused only once its
   *  plane lies below everything read by the current band, and its content is written back to view at that
   *  moment. lowest is the most negative offset (<= 0) at which g reads view along axis 0; reads above the
   *  band hit untouched planes and need no buffering. Reads past the local block land on ghosts, which are
   *  never written.
   *
   *  Extra memory: B - lowest planes for a band of B planes. Launches: one per band, plus one flush.
   *
   *  Unit test: ctest -R test-sweepassignment
   **/
  template <typename View, typename R, size_t NDim>
  void sweepAssign(const View &view, const R &g, const LayoutStruct<NDim> &layout, int lowest, SweepBuffer &buffer)
  {
    using T = typename View::value_type;
    using device::Idx;
    using SweepAssignDetail::at;
    const device::IdxArray<NDim> sizes = layout.getSizesInMemory();
    const Idx nGhosts = layout.getNGhosts();
    const Idx nPlanes = sizes[0];
    Idx planeArea = 1;
    for (size_t d = 1; d < NDim; ++d)
      planeArea *= sizes[d];

    const Idx band = sweepBandWidth(planeArea, nPlanes, buffer.getBandWidth());
    const Idx nSlots = std::min<Idx>(band - lowest, nPlanes);

    device::IdxArray<NDim> ringSizes = sizes;
    ringSizes[0] = nSlots;
    T *ringData = buffer.get<T>(nSlots * planeArea);
    const auto ring = device::apply(
        [&](auto... n) { return device::memory::NDViewUnmanaged<T, NDim>(ringData, n...); }, ringSizes);

    const device::IdxArray<NDim> zero{};
    for (Idx first = 0; first < nPlanes; first += band) {
      device::IdxArray<NDim> stops = sizes;
      stops[0] = std::min(band, nPlanes - first);
      device::iteration::foreach (
          "SweepAssign", zero, stops, DEVICE_LAMBDA(const device::IdxArray<NDim> &i) {
            const Idx plane = first + i[0];
            device::IdxArray<NDim> slot = i, site;
            slot[0] = plane % nSlots;
            for (size_t d = 0; d < NDim; ++d)
              site[d] = i[d] + nGhosts;
            site[0] = plane + nGhosts;
            if (plane >= nSlots) {
              device::IdxArray<NDim> done = site;
              done[0] -= nSlots;
              at(view, done) = at(ring, slot);
            }
            at(ring, slot) = device::apply([&](auto... j) { return DoEval::eval(g, j...); }, site);
          });
    }

    // Write back the planes still in the ring.
    device::IdxArray<NDim> stops = sizes;
    stops[0] = nSlots;
    device::iteration::foreach (
        "SweepAssignFlush", zero, stops, DEVICE_LAMBDA(const device::IdxArray<NDim> &i) {
          const Idx plane = nPlanes - nSlots + i[0];
          device::IdxArray<NDim> slot = i, site;
          slot[0] = plane % nSlots;
          for (size_t d = 0; d < NDim; ++d)
            site[d] = i[d] + nGhosts;
          site[0] = plane + nGhosts;
          at(view, site) = at(ring, slot);
        });
  }
} // namespace TempLat

#endif
