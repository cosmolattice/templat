
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026
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

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<2>> test2;
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<3>> test3;
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<4>> test4;
  TempLat::TDDContainer<TempLat::GhostUpdaterBatchTester<5>> test5;
} // namespace
