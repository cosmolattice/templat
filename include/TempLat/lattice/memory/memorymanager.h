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

    template <typename R = T>
    auto getNDSubView(const device::IdxArray<NDim> &localSizes,
                      const device::array<std::pair<device::Idx, device::Idx>, NDim> &slices) const
    {
      auto view = mBlock.template getNDView<R>(localSizes);
      auto subView = device::apply([&](const auto &...args) { return device::memory::subview(view, args...); }, slices);
      return subView;
    }
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
    device::memory::host_ptr<MemoryToolBox<NDim>> mToolBox;
    std::string mName;
    bool mAllocated;
    MemoryBlock<T, NDim> mBlock;

    MemoryLayoutState mLayoutState;
    GhostStateKeeper mGhostStateKeeper;
  };

} // namespace TempLat

#endif
