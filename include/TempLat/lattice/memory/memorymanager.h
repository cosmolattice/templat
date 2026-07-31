#ifndef TEMPLAT_LATTICE_MEMORY_MEMORYMANAGER_H
#define TEMPLAT_LATTICE_MEMORY_MEMORYMANAGER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/ghostcells/ghoststatekeeper.h"
#include "TempLat/lattice/memory/memoryblock.h"
#include "TempLat/lattice/memory/memorylayoutstate.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/parallel/device_memory.h"

#include <span>
#include <vector>
namespace TempLat
{
  MakeException(MemoryManagerAccessOutOfBounds);

  /** @brief A class which holds a single lattice in memory, and tracks and moves between various ghost states.
   *Templated for the memory type, typically float or double. All the confirm___Space() functions return a device::Idx,
   *which counts the amount of work done. Ignore that, but it is important for testing purposes. NOTE that the memory
   *access operators (T* and operator[]) do NOT verify if memory was allocated.
   *
   * Unit test: ctest -R test-memorymanager
   **/
  template <typename T, size_t NDim> class MemoryManager
  {
  public:
    // Put public methods here. These should change very little over time.
    MemoryManager(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, std::string name = "")
        : mToolBox(toolBox), mName(name), mAllocated(false)
    {
      // A null toolbox is not this class's error to report: ConfigView's constructor throws
      // FieldViewConfigMissingToolBox immediately after this one returns, exactly as it did
      // before the cold block moved here. Just leave the block default-constructed.
      if (mToolBox != nullptr) initConfigColdBlock();
    }

    /** @brief The configuration-space cold block, moved out of ConfigView.
     *
     * These three are pure functions of the toolbox's config-space layout, written once and
     * never again, and each used to be a member of ConfigView -- which is to say a member of
     * every leaf of every expression tree, deep-copied once per level of tree depth, because
     * expression nodes store their operands by value with the deduced type. 360 B of the
     * 456 B ConfigView was this. The manager is already one per field and already reference
     * counted, so holding them here costs nothing per copy.
     *
     * Returned by CONST REFERENCE, never by value: a 288 B layout copy on a path reachable
     * from expression evaluation is precisely the cost this removes.
     *
     * CONFIG-SCOPED IN THE NAME, DELIBERATELY. A Field's ConfigView and its FourierView share
     * ONE MemoryManager (field.h builds the Fourier view from *this, copying mManager), but
     * the two views hold DIFFERENT values for exactly these quantities:
     *
     *   layout         config-space layout           vs. Fourier-space layout
     *   memory sizes   real, ghost-padded            vs. Fourier, complex<T>-shaped
     *   local slicing  (nGhosts, nGhosts + local)    vs. never assigned -- no ghosts in k-space
     *
     * so one unscoped set here cannot serve both. The tempting symmetric cleanup -- having
     * FourierView::getFullNDHostView() read these too -- would hand it the config extents and
     * produce a host view with the wrong strides and silently wrong data; and the two slicing
     * types are not even interconvertible, so it need not fail to compile in a way that warns
     * you. FourierView keeps its own memorySizes. Do not unify these.
     */
    const LayoutStruct<NDim> &configLayout() const { return mConfigLayout; }
    const device::IdxArray<NDim> &configMemorySizes() const { return mConfigMemorySizes; }
    const device::array<std::pair<device::Idx, device::Idx>, NDim> &configLocalSlicing() const
    {
      return mConfigLocalSlicing;
    }

    device::Idx allocate()
    {
      if (mAllocated) return 0;
      const size_t size = mToolBox->mLayouts.getNecessaryMemoryAllocation();
      if (mToolBox->verbosity.allocation) sayMPI << "Allocating memory with " << size << " elements.\n";
      mAllocated = true;
      mBlock.allocate(size);
      return 1;
    }

    template <typename R = T> auto getNDView(const device::IdxArray<NDim> &localSizes) const
    {
      return mBlock.template getNDView<R>(localSizes);
    }
    template <typename R = T> auto getNDHostView(const device::IdxArray<NDim> &localSizes) const
    {
      return mBlock.template getNDHostView<R>(localSizes);
    }

    // getNDSubView (the device-side twin of getNDHostSubView below) is gone: it had no
    // callers anywhere in templat, CosmoInterface or MCInterface. Its host counterpart is
    // live -- ConfigView::getLocalNDHostView uses it.
    template <typename R = T>
    auto getNDHostSubView(const device::IdxArray<NDim> &localSizes,
                          const device::array<std::pair<device::Idx, device::Idx>, NDim> &slices) const
    {
      auto view = mBlock.template getNDHostView<R>(localSizes);
      auto subView = device::apply([&](const auto &...args) { return device::memory::subview(view, args...); }, slices);
      return subView;
    }

    void pushHostView() { mBlock.pushHostView(); }

    void deallocateHostView() { mBlock.deallocateHostView(); }

    template <typename R = T> auto getRawView() const { return mBlock.template getRawView<R>(); }

    template <typename R = T> auto getRawHostView() const { return mBlock.template getRawHostView<R>(); }

    device::Idx confirmConfigSpace()
    {
      device::Idx result = allocate();

      if (mToolBox->verbosity.spaceConfirmation)
        sayMPI << "Confirming that we are in configuration space. " << getName() << "\n";
      if (!mLayoutState.isConfigSpace()) {
        if (mLayoutState.isFourierSpace()) {
          if (mToolBox->verbosity.spaceConfirmation) sayMPI << "Need FFT C2R.\n";
          // say << "Setting fft library verbose.\n";
          // mToolBox->mFFTLibrary.setVerbose();
          // do an fft
          ++result;
          if (mToolBox->verbosity.fftPerformance) say << "FFT: " << mName << "(k) -> " << mName << "(x)\n";
          mToolBox->mFFTLibrary.c2r(mBlock);
          // normalize after FFT
          ++result;
          mToolBox->mFFTNormalization.c2r(mBlock, T{1});
          mLayoutState.setToFFTConfigSpace();
          if (mToolBox->verbosity.spaceConfirmation) sayMPI << "Performed FFT C2R.\n";
        }
        if (mLayoutState.isFFTConfigSpace()) {
          if (mToolBox->verbosity.spaceConfirmation)
            sayMPI << "Need ghost buster from fft config to plain config space.\n";
          ++result;
          mToolBox->mGhostBuster_toConfig(mBlock);
        }
        if (mToolBox->verbosity.spaceConfirmation)
          sayMPI << "Setting ghost state to stale, because of the FFT we performed. " << getName() << "\n";
        mGhostStateKeeper.setStale();
      }
      mLayoutState.setToConfigSpace();
      if (mToolBox->verbosity.spaceConfirmation) sayMPI << "We are in configuration space.\n";
      return result;
    }

    device::Idx confirmFFTConfigSpace()
    {
      device::Idx result = allocate();

      if (mToolBox->verbosity.spaceConfirmation)
        sayMPI << "Confirming that we are in FFT configuration space. " << getName() << "\n";
      if (!mLayoutState.isFFTConfigSpace()) {
        if (mLayoutState.isFourierSpace()) {
          if (mToolBox->verbosity.spaceConfirmation) sayMPI << "Need FFT C2R.\n";
          // say << "Setting fft library verbose.\n";
          // mToolBox->mFFTLibrary.setVerbose();
          // do an fft
          ++result;
          if (mToolBox->verbosity.fftPerformance) say << "FFT: " << mName << "(k) -> " << mName << "(x)\n";
          mToolBox->mFFTLibrary.c2r(mBlock);
          // normalize after FFT
          ++result;
          mToolBox->mFFTNormalization.c2r(mBlock, T{1});
          mLayoutState.setToFFTConfigSpace();
          if (mToolBox->verbosity.spaceConfirmation) sayMPI << "Performed FFT C2R.\n";
        } else if (mLayoutState.isConfigSpace()) {
          if (mToolBox->verbosity.spaceConfirmation)
            sayMPI << "Need ghost buster from plain config to fft config space.\n";
          ++result;
          mToolBox->mGhostBuster_toFFTConfig(mBlock);
          mLayoutState.setToFFTConfigSpace();
        }
      }
      mLayoutState.setToFFTConfigSpace();
      if (mToolBox->verbosity.spaceConfirmation) sayMPI << "We are in FFT configuration space.\n";
      return result;
    }

    device::Idx confirmFourierSpace()
    {
      device::Idx result = allocate();

      if (mToolBox->verbosity.spaceConfirmation)
        sayMPI << "Confirming that we are in Fourier space. " << getName() << "\n";
      if (!mLayoutState.isFourierSpace()) {
        if (mLayoutState.isConfigSpace()) {
          if (mToolBox->verbosity.spaceConfirmation)
            sayMPI << "Need ghost buster from plain config to fft config space.\n";
          ++result;
          mToolBox->mGhostBuster_toFFTConfig(mBlock);
          mLayoutState.setToFFTConfigSpace();
        }
        if (mLayoutState.isFFTConfigSpace()) {
          if (mToolBox->verbosity.spaceConfirmation) sayMPI << "Need FFT R2C.\n";
          ++result;
          if (mToolBox->verbosity.fftPerformance) say << "FFT: " << mName << "(x) -> " << mName << "(k)\n";
          // do an fft
          mToolBox->mFFTLibrary.r2c(mBlock);
          // normalize after FFT
          ++result;
          mToolBox->mFFTNormalization.r2c(mBlock, T{1});
          if (mToolBox->verbosity.spaceConfirmation) sayMPI << "Performed FFT R2C.\n";
        }
        mBlock.flagHostMirrorOutdated();
      }
      mLayoutState.setToFourierSpace();
      if (mToolBox->verbosity.spaceConfirmation) sayMPI << "We are in Fourier space.\n";
      return result;
    }

    device::Idx confirmGhostsUpToDate()
    {
      device::Idx result = confirmConfigSpace();

      if (mToolBox->verbosity.ghostConfirmationSteps)
        sayMPI << "Confirming that ghost cells are up to date. " << getName() << "\n" << mGhostStateKeeper << "\n";

      if (mGhostStateKeeper.isStale()) {
        if (mToolBox->verbosity.ghostConfirmationSteps) sayMPI << "Need to update ghost cells.\n";
        ++result;
        mToolBox->mGhostUpdater.update(mBlock);
        mGhostStateKeeper.setUpToDate();
      }
      if (mToolBox->verbosity.ghostConfirmationSteps)
        sayMPI << "Ghost cells are up to date.\n" << mGhostStateKeeper << "\n";
      mBlock.flagHostMirrorOutdated();
      return result;
    }

    void flagHostMirrorOutdated() { mBlock.flagHostMirrorOutdated(); }

    void updateGhosts()
    {
      mToolBox->mGhostUpdater.update(mBlock);
      mBlock.flagHostMirrorOutdated();
    }

    /** @brief Coalesced confirmGhostsUpToDate over the managers of a multi-component field. Only the
     *  components whose ghosts are stale are exchanged, but they travel as one message per
     *  dimension/direction instead of one per component. Static so it may reach the private state of the
     *  sibling managers; they must all share one MemoryToolBox (one GhostUpdater / buffer set). */
    static device::Idx confirmGhostsUpToDateBatch(std::span<MemoryManager<T, NDim> *const> mgrs)
    {
      device::Idx result = 0;
      if (mgrs.empty()) return result;
      auto toolBox = mgrs[0]->mToolBox;
      std::vector<MemoryBlock<T, NDim> *> stale;
      stale.reserve(mgrs.size());
      for (auto *m : mgrs) {
        if (m->mToolBox.get() != toolBox.get())
          throw MemoryManagerAccessOutOfBounds("confirmGhostsUpToDateBatch: managers do not share a MemoryToolBox.");
        result += m->confirmConfigSpace();
        if (m->mGhostStateKeeper.isStale()) stale.push_back(&m->mBlock);
      }
      if (!stale.empty()) {
        ++result;
        toolBox->mGhostUpdater.updateBatch(std::span<MemoryBlock<T, NDim> *const>(stale.data(), stale.size()));
        for (auto *m : mgrs)
          m->mGhostStateKeeper.setUpToDate();
      }
      for (auto *m : mgrs)
        m->mBlock.flagHostMirrorOutdated();
      return result;
    }

    /** @brief Coalesced forced updateGhosts over the managers of a multi-component field. Matches the
     *  per-component updateGhosts() semantics (unconditional, no state-keeper reset) but sends all
     *  components in one message per dimension/direction. */
    static void updateGhostsBatch(std::span<MemoryManager<T, NDim> *const> mgrs)
    {
      if (mgrs.empty()) return;
      auto toolBox = mgrs[0]->mToolBox;
      std::vector<MemoryBlock<T, NDim> *> blocks;
      blocks.reserve(mgrs.size());
      for (auto *m : mgrs) {
        if (m->mToolBox.get() != toolBox.get())
          throw MemoryManagerAccessOutOfBounds("updateGhostsBatch: managers do not share a MemoryToolBox.");
        blocks.push_back(&m->mBlock);
      }
      toolBox->mGhostUpdater.updateBatch(std::span<MemoryBlock<T, NDim> *const>(blocks.data(), blocks.size()));
      for (auto *m : mgrs)
        m->mBlock.flagHostMirrorOutdated();
    }

    /** @brief this is the only state the one may need to set from the outside: if a field is updated in the integrator.
     */
    void setGhostsAreStale()
    {
      if (mToolBox->verbosity.ghostConfirmationSteps)
        sayMPI << "Noting that ghost cells are no longer up to date. " << getName() << "\n";
      mGhostStateKeeper.setStale();
    }

    bool areGhostsStale() const { return mGhostStateKeeper.isStale(); }

    /** @brief Check the current state. */
    bool isConfigSpace() const { return mLayoutState.isConfigSpace(); }

    /** @brief Check the current state. */
    bool isFourierSpace() const { return mLayoutState.isFourierSpace(); }

    std::string getName() const { return mName; }

    void setName(std::string newName) { mName = newName; }

    friend std::ostream &operator<<(std::ostream &ostream, const MemoryManager &mMan)
    {
      ostream << "Memory manager -\n " << mMan.mToolBox->mGroup << "\n"
              << mMan.mToolBox->mLayouts << "\n\n"
              << mMan.mLayoutState << "\n\n"
              << mMan.mGhostStateKeeper;
      return ostream;
    }

    friend bool operator==(const MemoryManager &a, const MemoryManager &b) { return a.mBlock == b.mBlock; }

    size_t bytes() const { return mBlock.size() * sizeof(T); }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    /** @brief Derive the config-space cold block from the toolbox layouts.
     *
     * Verbatim what ConfigView's constructor used to do. It is purely derived from the
     * layout -- no allocation is involved -- so running it this early, before the first
     * confirmConfigSpace(), is safe. What is NOT safe to move earlier is ConfigView's
     * mView: that caches mData.data() from the MemoryBlock and so must stay after the
     * confirmConfigSpace() that allocates it.
     */
    void initConfigColdBlock()
    {
      mConfigLayout = mToolBox->mLayouts.getConfigSpaceLayout();

      const auto localSizes = mConfigLayout.getLocalSizes();
      const size_t nGhosts = mConfigLayout.getNGhosts();

      mConfigMemorySizes = mConfigLayout.getSizesInMemory();
      for (size_t d = 0; d < NDim; ++d) {
        mConfigMemorySizes[d] += nGhosts + nGhosts; // add padding to the local sizes
        mConfigLocalSlicing[d] = std::make_pair(nGhosts, nGhosts + localSizes[d]);
      }
    }

    device::memory::host_ptr<MemoryToolBox<NDim>> mToolBox;
    std::string mName;
    bool mAllocated;
    MemoryBlock<T, NDim> mBlock;

    MemoryLayoutState mLayoutState;
    GhostStateKeeper mGhostStateKeeper;

    // The config-scoped cold block -- see configLayout() above for why the scoping is
    // load-bearing and must not be shared with FourierView.
    LayoutStruct<NDim> mConfigLayout;
    device::IdxArray<NDim> mConfigMemorySizes;
    device::array<std::pair<device::Idx, device::Idx>, NDim> mConfigLocalSlicing;
  };

} // namespace TempLat

#endif
