
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
//
// Verifies the coalesced multi-block ghost exchange (GhostUpdater::updateBatch, used by the
// component-coalesced updateGhosts of ComplexField / SU2Field / SymTracelessField). Each of C blocks
// is filled with its global coordinate AND its component index; after one updateBatch every padded
// cell (corners included) must equal the periodic-wrapped coordinate and still carry its OWN component
// label — a mis-ordered coalesced pack would either corrupt the coordinate or swap component labels.

#include "TempLat/lattice/ghostcells/ghostupdater.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/lattice/memory/triplestatelayouts.h"

#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/util/log/saycomplete.h"

#include <cmath>
#include <sstream>
#include <vector>
#include <memory>
#include <span>

namespace TempLat
{
  namespace TestScratch
  {
    template <size_t NDim> struct labelledMPITypeHolder {
      MPI_Datatype dType;
      labelledMPITypeHolder()
      {
#ifdef HAVE_MPI
        MPI_Type_contiguous(NDim + 1, TempLat::MPITypeSelect<device::Idx>(), &dType);
        MPI_Type_commit(&dType);
#endif
      }
      ~labelledMPITypeHolder()
      {
#ifdef HAVE_MPI
        int didFinalize = 0;
        MPI_Finalized(&didFinalize);
        if (!didFinalize) MPI_Type_free(&dType);
#endif
      }
    };

    template <size_t NDim> struct labelled {
      device::IdxArray<NDim> data;
      device::Idx component;

      static MPI_Datatype getMPIType()
      {
        static labelledMPITypeHolder<NDim> holder;
        return holder.dType;
      }
    };

    template <size_t NDim>
    void labelled_initialize(MemoryBlock<labelled<NDim>, NDim> &block, const LayoutStruct<NDim> layout,
                             device::Idx component)
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
      auto subView = device::apply([&](const auto &...args) { return device::memory::subview(view, args...); }, slices);

      auto functor = DEVICE_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        labelled<NDim> val;
        for (size_t k = 0; k < NDim; ++k)
          val.data[k] = idx[k] + layout.getLocalStarts()[k] + 1;
        val.component = component;
        device::IdxArray<NDim> inc;
        for (size_t k = 0; k < NDim; ++k)
          inc[k] = k;
        device::apply([&](const auto &...args) { subView(idx[args]...) = val; }, inc);
      };

      device::IdxArray<NDim> it_stop{};
      for (size_t k = 0; k < NDim; ++k)
        it_stop[k] = localSizes[k];
      device::iteration::foreach ("LabelledInit", {}, it_stop, functor);
      device::iteration::fence();
    }

    /** @brief Fills the owned region with a value that is distinct per (site, component). A constant
     *  field would not distinguish Neumann's outward mirror from a plain copy, and identical values
     *  across components would hide cross-component corruption in the coalesced pack. */
    template <size_t NDim>
    void double_initialize(MemoryBlock<double, NDim> &block, const LayoutStruct<NDim> layout,
                           device::Idx component)
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
      auto subView = device::apply([&](const auto &...args) { return device::memory::subview(view, args...); }, slices);

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
      device::iteration::foreach ("DoubleInit", {}, it_stop, functor);
      device::iteration::fence();
    }

    /** @brief The coalesced batch exchange must be indistinguishable from exchanging each component
     *  on its own. Runs both paths over identically-initialized blocks and compares every padded cell,
     *  ghosts and corners included.
     *
     *  Regression guard: updateBatch() originally took no BCSpec, so the batch path silently applied
     *  Periodic no matter what BC the field carried. Every non-Periodic case below diverged from the
     *  single-block path — which is the path a scalar Field takes, and therefore the one all the other
     *  BC tests in this tree cover. */
    template <size_t nd>
    bool test_batch_matches_single(const device::Idx nGrid, const size_t nGhost, const size_t C,
                                   size_t bcDim, BCType bc)
    {
      auto toolBox = MemoryToolBox<nd>::makeShared(nGrid, nGhost);
      toolBox->unsetVerbose();

      MPICartesianGroup mGroup(FFTMPIDomainSplit<nd>::makeMPIGroup(toolBox->mNGridPointsVec));
      FFTLibrarySelector<nd> fftlib(mGroup, toolBox->mNGridPointsVec);
      TripleStateLayouts fullLayout(fftlib.getLayout(), nGhost);
      GhostUpdater<nd> ghostUpdater(mGroup, fullLayout.getConfigSpaceLayout());
      const auto config_layout = fullLayout.getConfigSpaceLayout();

      const auto localSizes = config_layout.getLocalSizes();
      size_t total_size = 1;
      for (size_t i = 0; i < nd; ++i)
        total_size *= localSizes[i] + 2 * nGhost;

      BCSpec<nd> spec = allPeriodic<nd>();
      spec[bcDim] = bc;

      std::vector<std::unique_ptr<MemoryBlock<double, nd>>> batched, single;
      std::vector<MemoryBlock<double, nd> *> ptrs;
      batched.reserve(C);
      single.reserve(C);
      ptrs.reserve(C);
      for (size_t c = 0; c < C; ++c) {
        batched.push_back(std::make_unique<MemoryBlock<double, nd>>(total_size));
        single.push_back(std::make_unique<MemoryBlock<double, nd>>(total_size));
        double_initialize<nd>(*batched[c], config_layout, static_cast<device::Idx>(c));
        double_initialize<nd>(*single[c], config_layout, static_cast<device::Idx>(c));
        ptrs.push_back(batched[c].get());
      }

      ghostUpdater.updateBatch(std::span<MemoryBlock<double, nd> *const>(ptrs.data(), ptrs.size()), spec);
      for (size_t c = 0; c < C; ++c)
        ghostUpdater.update(*single[c], spec);

      bool ok = true;
      for (size_t c = 0; c < C; ++c) {
        batched[c]->flagHostMirrorOutdated();
        single[c]->flagHostMirrorOutdated();
        auto vb = batched[c]->getRawHostView();
        auto vs = single[c]->getRawHostView();
        for (size_t i = 0; i < total_size; ++i) {
          if (std::abs(vb(i) - vs(i)) > 1e-12) {
            ok = false;
            std::stringstream ss;
            ss << "batch/single mismatch: nd=" << nd << " bc=" << static_cast<int>(bc) << " bcDim=" << bcDim
               << " comp=" << c << " flat=" << i << " batch=" << vb(i) << " single=" << vs(i) << "\n";
            sayMPI << ss.str();
            break; // one report per component is enough to localize the failure
          }
        }
      }
      return ok;
    }

    template <size_t nd> bool test_ghost_updater_batch(const device::Idx nGrid, const size_t nGhost, const size_t C)
    {
      auto toolBox = MemoryToolBox<nd>::makeShared(nGrid, nGhost);
      toolBox->unsetVerbose();

      MPICartesianGroup mGroup(FFTMPIDomainSplit<nd>::makeMPIGroup(toolBox->mNGridPointsVec));
      FFTLibrarySelector<nd> fftlib(mGroup, toolBox->mNGridPointsVec);
      TripleStateLayouts fullLayout(fftlib.getLayout(), nGhost);
      GhostUpdater<nd> ghostUpdater(mGroup, fullLayout.getConfigSpaceLayout());
      const auto config_layout = fullLayout.getConfigSpaceLayout();

      const auto localSizes = config_layout.getLocalSizes();
      device::IdxArray<nd> full{};
      size_t total_size = 1;
      for (size_t i = 0; i < nd; ++i) {
        full[i] = localSizes[i] + 2 * nGhost;
        total_size *= full[i];
      }

      // MemoryBlock is non-movable, so hold the C blocks via unique_ptr.
      std::vector<std::unique_ptr<MemoryBlock<labelled<nd>, nd>>> blocks;
      std::vector<MemoryBlock<labelled<nd>, nd> *> ptrs;
      blocks.reserve(C);
      ptrs.reserve(C);
      for (size_t c = 0; c < C; ++c) {
        blocks.push_back(std::make_unique<MemoryBlock<labelled<nd>, nd>>(total_size));
        labelled_initialize<nd>(*blocks[c], config_layout, static_cast<device::Idx>(c));
        ptrs.push_back(blocks[c].get());
      }

      ghostUpdater.updateBatch(std::span<MemoryBlock<labelled<nd>, nd> *const>(ptrs.data(), ptrs.size()));

      bool all_correct = true;
      for (size_t c = 0; c < C; ++c) {
        blocks[c]->flagHostMirrorOutdated();
        auto view = blocks[c]->getRawHostView();
        device::IdxArray<nd> cIdx{};
        for (size_t i = 0; i < total_size; ++i) {
          size_t remainder = i;
          for (size_t j = 0; j < nd; ++j) {
            cIdx[nd - 1 - j] = remainder % full[nd - 1 - j];
            remainder /= full[nd - 1 - j];
            cIdx[nd - 1 - j] += config_layout.getLocalStarts()[nd - 1 - j] + 1 - (device::Idx)nGhost;
          }
          auto should = cIdx;
          for (size_t d = 0; d < nd; ++d) {
            if (cIdx[d] < 1) should[d] += nGrid;
            if (cIdx[d] > (device::Idx)nGrid) should[d] -= nGrid;
          }
          for (size_t d = 0; d < nd; ++d)
            all_correct &= (size_t)view(i).data[d] == (size_t)should[d];
          all_correct &= (size_t)view(i).component == c; // component label preserved -> no cross-talk
        }
      }
      return all_correct;
    }
  } // namespace TestScratch

  template <size_t NDim> struct GhostUpdaterBatchTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void GhostUpdaterBatchTester<NDim>::Test(TDDAssertion &tdd)
  {
    static_assert(NDim > 1, "GhostUpdater batch test only makes sense in 2 or more dimensions.");
    // Component counts mirror the real fields: 2 (complex), 4 (SU2), 5 (symmetric traceless).
    if constexpr (NDim < 6) {
      tdd.verify(TestScratch::test_ghost_updater_batch<NDim>(16, 1, 2));
      tdd.verify(TestScratch::test_ghost_updater_batch<NDim>(16, 1, 4));
      tdd.verify(TestScratch::test_ghost_updater_batch<NDim>(16, 2, 5));
    }
    if constexpr (NDim < 4) {
      tdd.verify(TestScratch::test_ghost_updater_batch<NDim>(32, 1, 4));
      tdd.verify(TestScratch::test_ghost_updater_batch<NDim>(32, 2, 2));
    }
  }

  template <size_t NDim> struct GhostUpdaterBatchBCTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void GhostUpdaterBatchBCTester<NDim>::Test(TDDAssertion &tdd)
  {
    static_assert(NDim > 1, "GhostUpdater batch test only makes sense in 2 or more dimensions.");
    // Sweep every BC over every dimension: dim 0 is the MPI-split one under the usual decomposition,
    // the rest exercise the non-split local path inside pUpdateBatch.
    for (auto bc : {BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet, BCType::Neumann}) {
      for (size_t d = 0; d < NDim; ++d) {
        tdd.verify(TestScratch::test_batch_matches_single<NDim>(16, 1, 2, d, bc));
        tdd.verify(TestScratch::test_batch_matches_single<NDim>(16, 1, 4, d, bc));
      }
      // Deeper halo: the BC post-step loops depth 1..nGhost, so nGhost > 1 is a distinct path.
      tdd.verify(TestScratch::test_batch_matches_single<NDim>(16, 2, 5, 0, bc));
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<2>> test2;
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<3>> test3;
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<4>> test4;
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<5>> test5;

  TempLat::TDDContainer<TempLat::GhostUpdaterBatchBCTester<2>> bcTest2;
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchBCTester<3>> bcTest3;
} // namespace
