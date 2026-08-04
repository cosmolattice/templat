
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026
//
// GhostUpdater::updateBatch's PER-COMPONENT BCSpec overload, one level below the MemoryManager
// entry points that percomponentbc.cpp drives.
//
//   updateBatch(blocks, std::span<const BCSpec<NDim>>)   -- bcSpecs[c] belongs to blocks[c]
//
// The contract, from the top down:
//   * the coalesced exchange is BC-agnostic: it lands the global-wrap value in every ghost slab,
//     and the BC is a per-block post-step, so components may legitimately disagree (this is what
//     C-star is -- Phi(x+L) = Phi*(x) flips c1 and c3 of an SU(2) element but not c0 and c2);
//   * the fused LOCAL kernel writes one BCType per launch, so the batch is partitioned by
//     BCType-in-this-dimension and launched once per distinct type present -- a uniform batch keeps
//     its single launch, a C-star batch costs two;
//   * the older single-BCSpec overload broadcasts into this one, so every pre-existing caller is
//     source-compatible and a uniform batch is unchanged.
//
// Launch counts are not observable from here, so the partition is covered behaviourally: whatever
// the grouping, the result must equal updating each block on its own with its own BCSpec, cell for
// cell. That comparison also covers the broadcast (a uniform span and the single-BCSpec overload
// must produce identical bytes).

#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/ghostcells/ghostupdater.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/lattice/memory/triplestatelayouts.h"
#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/util/log/saycomplete.h"
#include "TempLat/util/tdd/tdd.h"

#include <memory>
#include <span>
#include <sstream>
#include <vector>

namespace TempLat
{
  namespace PerComponentBCScratch
  {
    /** @brief Distinct per (site, component): a constant fill would not tell Neumann's outward
     *  mirror from a plain copy, and identical values across components would hide a coalesced pack
     *  that swapped two components' slabs. Varies in every dimension so a corner ghost taken from
     *  the wrong dimension is visible too. */
    template <size_t NDim>
    void labelledInit(MemoryBlock<double, NDim> &block, const LayoutStruct<NDim> layout, device::Idx component)
    {
      const auto localSizes = layout.getLocalSizes();
      const device::Idx nGhost = layout.getNGhosts();

      device::IdxArray<NDim> viewSizes;
      for (size_t k = 0; k < NDim; ++k)
        viewSizes[k] = localSizes[k] + 2 * nGhost;
      auto view = block.getNDView(viewSizes);

      device::array<std::pair<device::Idx, device::Idx>, NDim> slices{};
      for (size_t k = 0; k < NDim; ++k)
        slices[k] = std::make_pair(nGhost, nGhost + localSizes[k]);
      auto subView =
          device::apply([&](const auto &...args) { return device::memory::subview(view, args...); }, slices);

      auto functor = DEVICE_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        double val = 1.0e9 * (double)(component + 1);
        double mult = 1.0;
        for (size_t k = 0; k < NDim; ++k) {
          val += (double)(idx[k] + layout.getLocalStarts()[k] + 1) * mult;
          mult *= 37.0;
        }
        device::IdxArray<NDim> inc;
        for (size_t k = 0; k < NDim; ++k)
          inc[k] = k;
        device::apply([&](const auto &...args) { subView(idx[args]...) = val; }, inc);
      };

      device::IdxArray<NDim> it_stop{};
      for (size_t k = 0; k < NDim; ++k)
        it_stop[k] = localSizes[k];
      device::iteration::foreach ("PerComponentBCInit", {}, it_stop, functor);
      device::iteration::fence();
    }

    /** @brief Owns a GhostUpdater plus N identically-initialized blocks, so several update paths can
     *  be run over the same starting data and compared. */
    template <size_t NDim> struct Harness {
      device::memory::host_ptr<MemoryToolBox<NDim>> toolBox;
      MPICartesianGroup group;
      FFTLibrarySelector<NDim> fftlib;
      TripleStateLayouts<NDim> layouts;
      GhostUpdater<NDim> ghostUpdater;
      LayoutStruct<NDim> configLayout;
      size_t totalSize;

      Harness(device::Idx nGrid, size_t nGhost)
          : toolBox(MemoryToolBox<NDim>::makeShared(nGrid, nGhost)),
            group(FFTMPIDomainSplit<NDim>::makeMPIGroup(toolBox->mNGridPointsVec)),
            fftlib(group, toolBox->mNGridPointsVec), layouts(fftlib.getLayout(), nGhost),
            ghostUpdater(group, layouts.getConfigSpaceLayout()), configLayout(layouts.getConfigSpaceLayout())
      {
        toolBox->unsetVerbose();
        const auto localSizes = configLayout.getLocalSizes();
        totalSize = 1;
        for (size_t i = 0; i < NDim; ++i)
          totalSize *= localSizes[i] + 2 * nGhost;
      }

      std::vector<std::unique_ptr<MemoryBlock<double, NDim>>> makeBlocks(size_t C) const
      {
        std::vector<std::unique_ptr<MemoryBlock<double, NDim>>> blocks;
        blocks.reserve(C);
        for (size_t c = 0; c < C; ++c) {
          blocks.push_back(std::make_unique<MemoryBlock<double, NDim>>(totalSize));
          labelledInit<NDim>(*blocks.back(), configLayout, static_cast<device::Idx>(c));
        }
        return blocks;
      }

      bool identical(std::vector<std::unique_ptr<MemoryBlock<double, NDim>>> &a,
                     std::vector<std::unique_ptr<MemoryBlock<double, NDim>>> &b, const char *what) const
      {
        bool ok = true;
        for (size_t c = 0; c < a.size(); ++c) {
          a[c]->flagHostMirrorOutdated();
          b[c]->flagHostMirrorOutdated();
          auto va = a[c]->getRawHostView();
          auto vb = b[c]->getRawHostView();
          for (size_t i = 0; i < totalSize; ++i) {
            if (va(i) != vb(i)) {
              ok = false;
              std::stringstream ss;
              ss << "[" << what << "] mismatch: NDim=" << NDim << " comp=" << c << " flat=" << i << " lhs=" << va(i)
                 << " rhs=" << vb(i) << "\n";
              sayMPI << ss.str();
              break;
            }
          }
        }
        return ok;
      }
    };

    template <size_t NDim>
    std::vector<MemoryBlock<double, NDim> *> raw(std::vector<std::unique_ptr<MemoryBlock<double, NDim>>> &blocks)
    {
      std::vector<MemoryBlock<double, NDim> *> ptrs;
      ptrs.reserve(blocks.size());
      for (auto &b : blocks)
        ptrs.push_back(b.get());
      return ptrs;
    }

    /** @brief The per-component span overload must agree, cell for cell, with updating each block on
     *  its own under its own BCSpec. Holds whatever the internal partition by BCType does. */
    template <size_t NDim>
    bool batchSpanMatchesSingle(device::Idx nGrid, size_t nGhost, const std::vector<BCSpec<NDim>> &specs)
    {
      Harness<NDim> h(nGrid, nGhost);
      auto batched = h.makeBlocks(specs.size());
      auto single = h.makeBlocks(specs.size());
      auto ptrs = raw<NDim>(batched);

      h.ghostUpdater.updateBatch(std::span<MemoryBlock<double, NDim> *const>(ptrs.data(), ptrs.size()),
                                 std::span<const BCSpec<NDim>>(specs.data(), specs.size()));
      for (size_t c = 0; c < specs.size(); ++c)
        h.ghostUpdater.update(*single[c], specs[c]);

      return h.identical(batched, single, "span-vs-single");
    }

    /** @brief Source compatibility and the preserved fast path in one: for a uniform batch the old
     *  single-BCSpec overload, the span overload, and the per-component single-block path must all
     *  produce identical bytes. The single-BCSpec overload is the one every pre-existing caller
     *  uses; it must still compile and still broadcast. */
    template <size_t NDim>
    bool uniformBroadcastIsUnchanged(device::Idx nGrid, size_t nGhost, size_t C, BCSpec<NDim> spec)
    {
      Harness<NDim> h(nGrid, nGhost);
      auto viaOldOverload = h.makeBlocks(C);
      auto viaSpan = h.makeBlocks(C);
      auto viaSingle = h.makeBlocks(C);

      auto pOld = raw<NDim>(viaOldOverload);
      auto pSpan = raw<NDim>(viaSpan);
      const std::vector<BCSpec<NDim>> specs(C, spec);

      h.ghostUpdater.updateBatch(std::span<MemoryBlock<double, NDim> *const>(pOld.data(), pOld.size()), spec);
      h.ghostUpdater.updateBatch(std::span<MemoryBlock<double, NDim> *const>(pSpan.data(), pSpan.size()),
                                 std::span<const BCSpec<NDim>>(specs.data(), specs.size()));
      for (size_t c = 0; c < C; ++c)
        h.ghostUpdater.update(*viaSingle[c], spec);

      bool ok = h.identical(viaOldOverload, viaSpan, "broadcast-vs-span");
      ok &= h.identical(viaOldOverload, viaSingle, "broadcast-vs-single");
      return ok;
    }
  } // namespace PerComponentBCScratch

  /** @brief updateBatch(blocks, span<BCSpec>) with components that disagree on their BC.
   *
   * Patterns: the SU(2)/doublet C-star pattern (P, AP, P, AP) and the complex scalar's (P, AP) --
   * two partition groups; all four BCTypes at once -- four groups, the widest partition the fused
   * kernel can be asked for; specs differing in different dimensions -- the partition is decided per
   * dimension, independently; and a batch wider than the fused kernel's capture budget, which falls
   * back to the per-component path and must index the BCSpecs there too.
   */
  template <size_t NDim> struct GhostUpdaterPerComponentBCTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void GhostUpdaterPerComponentBCTester<NDim>::Test(TDDAssertion &tdd)
  {
    static_assert(NDim > 1, "The batched ghost path only makes sense in 2 or more dimensions.");
    using namespace PerComponentBCScratch;

    auto withBC = [](size_t d, BCType bc) {
      BCSpec<NDim> s = allPeriodic<NDim>();
      s[d] = bc;
      return s;
    };

    for (size_t nGhost : {1u, 2u}) {
      for (size_t d = 0; d < NDim; ++d) {
        const std::vector<BCSpec<NDim>> cstar = {withBC(d, BCType::Periodic), withBC(d, BCType::Antiperiodic),
                                                 withBC(d, BCType::Periodic), withBC(d, BCType::Antiperiodic)};
        const std::vector<BCSpec<NDim>> complexCStar = {withBC(d, BCType::Periodic), withBC(d, BCType::Antiperiodic)};
        const std::vector<BCSpec<NDim>> allFour = {withBC(d, BCType::Periodic), withBC(d, BCType::Antiperiodic),
                                                   withBC(d, BCType::Dirichlet), withBC(d, BCType::Neumann)};
        tdd.verify(batchSpanMatchesSingle<NDim>(16, nGhost, cstar));
        tdd.verify(batchSpanMatchesSingle<NDim>(16, nGhost, complexCStar));
        tdd.verify(batchSpanMatchesSingle<NDim>(16, nGhost, allFour));
      }

      // Different dimensions disagree differently.
      std::vector<BCSpec<NDim>> crossDim(5, allPeriodic<NDim>());
      crossDim[1][0] = BCType::Antiperiodic;
      crossDim[2][NDim - 1] = BCType::Dirichlet;
      crossDim[3][0] = BCType::Antiperiodic;
      crossDim[3][NDim - 1] = BCType::Neumann;
      crossDim[4][0] = BCType::Neumann;
      tdd.verify(batchSpanMatchesSingle<NDim>(16, nGhost, crossDim));

      // Wider than the fused kernel's capture budget (cMaxLocalGhostBatch): per-component fallback.
      std::vector<BCSpec<NDim>> wide;
      for (size_t c = 0; c < 9; ++c)
        wide.push_back(withBC(0, (c % 2) ? BCType::Antiperiodic : BCType::Periodic));
      tdd.verify(batchSpanMatchesSingle<NDim>(16, nGhost, wide));
    }
  }

  /** @brief A uniform batch must be untouched by the per-component plumbing, and the single-BCSpec
   * overload every pre-existing caller uses must still broadcast into it. */
  template <size_t NDim> struct GhostUpdaterUniformBatchUnchangedTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void GhostUpdaterUniformBatchUnchangedTester<NDim>::Test(TDDAssertion &tdd)
  {
    static_assert(NDim > 1, "The batched ghost path only makes sense in 2 or more dimensions.");
    using namespace PerComponentBCScratch;

    for (auto bc : {BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet, BCType::Neumann}) {
      for (size_t d = 0; d < NDim; ++d) {
        BCSpec<NDim> spec = allPeriodic<NDim>();
        spec[d] = bc;
        tdd.verify(uniformBroadcastIsUnchanged<NDim>(16, 1, 4, spec));
      }
      BCSpec<NDim> spec = allPeriodic<NDim>();
      spec[0] = bc;
      tdd.verify(uniformBroadcastIsUnchanged<NDim>(16, 2, 5, spec));
    }
  }

  /** @brief One BCSpec per block, no more and no fewer: a span of the wrong length is a caller bug
   * that must not be silently padded or truncated. */
  template <size_t NDim> struct GhostUpdaterBCSpanSizeTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void GhostUpdaterBCSpanSizeTester<NDim>::Test(TDDAssertion &tdd)
  {
    using namespace PerComponentBCScratch;
    Harness<NDim> h(16, 1);
    auto blocks = h.makeBlocks(4);
    auto ptrs = raw<NDim>(blocks);
    const std::vector<BCSpec<NDim>> tooFew(3, allPeriodic<NDim>());
    const std::vector<BCSpec<NDim>> tooMany(5, allPeriodic<NDim>());

    tdd.verify(Throws<GhostUpdaterException>([&] {
      h.ghostUpdater.updateBatch(std::span<MemoryBlock<double, NDim> *const>(ptrs.data(), ptrs.size()),
                                 std::span<const BCSpec<NDim>>(tooFew.data(), tooFew.size()));
    }));
    tdd.verify(Throws<GhostUpdaterException>([&] {
      h.ghostUpdater.updateBatch(std::span<MemoryBlock<double, NDim> *const>(ptrs.data(), ptrs.size()),
                                 std::span<const BCSpec<NDim>>(tooMany.data(), tooMany.size()));
    }));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GhostUpdaterPerComponentBCTester<2>> guPerComponentBC2;
  TempLat::TDDContainer<TempLat::GhostUpdaterPerComponentBCTester<3>> guPerComponentBC3;

  TempLat::TDDContainer<TempLat::GhostUpdaterUniformBatchUnchangedTester<2>> guUniform2;
  TempLat::TDDContainer<TempLat::GhostUpdaterUniformBatchUnchangedTester<3>> guUniform3;

  TempLat::TDDContainer<TempLat::GhostUpdaterBCSpanSizeTester<3>> guSpanSize3;
} // namespace
