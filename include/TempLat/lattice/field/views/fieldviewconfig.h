#ifndef TEMPLAT_LATTICE_FIELD_VIEWS_FIELDVIEWCONFIG_H
#define TEMPLAT_LATTICE_FIELD_VIEWS_FIELDVIEWCONFIG_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/lattice/field/abstractfield.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstruct.h"

#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"

#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/parallel/device_iteration.h"

namespace TempLat
{
  MakeException(FieldViewConfigWrongSpaceConfirmation);
  MakeException(FieldViewConfigMissingToolBox);

  /** @brief A view on the field which, when interacted with, assures every time again that things are in
   *   configuration space, and possibly the ghost cells are updated when needed.
   *   The final Field class defaults to config space, which means it inherits from this class.
   *
   * Unit test: ctest -R test-fieldviewconfig
   **/
  template <typename T, size_t _NDim> class ConfigView : public AbstractField<T, _NDim>
  {
  public:
    // Put public methods here. These should change very little over time.
    static constexpr size_t NDim = _NDim;

    using AbstractField<T, NDim>::mManager;
    using AbstractField<T, NDim>::mToolBox;

    ConfigView(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, LatticeParameters<T> pLatPar)
        : AbstractField<T, NDim>(name, toolBox, pLatPar), mDisableFFTBlocking(false)
    {
      if (toolBox == nullptr)
        throw FieldViewConfigMissingToolBox("A FieldViewConfig must be constructed with a valid MemoryToolBox.");

      // The layout, the padded memory sizes and the local slicing now live in the
      // MemoryManager, which derived them from the same toolbox in its own constructor a
      // moment ago -- see MemoryManager::configLayout().

      mManager->setGhostsAreStale();
      mManager->confirmConfigSpace(); // allocation happens here

      // Order is load-bearing: mView caches mData.data() out of the MemoryBlock, and
      // confirmConfigSpace() above is what allocates it. Do not hoist this.
      mView = mManager->getNDView(mManager->configMemorySizes());
    }

    auto getView() const { return mView; }

    template <typename R> void assign(R &&g)
    {
      onBeforeAssignment(g);

      PreGet::apply(g);

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply([&](auto &&...args) { mView(args...) = DoEval::eval(g, args...); }, idx);
      };
      device::iteration::foreach ("ConfigViewAssign", mManager->configLayout(), functor);

      PostGet::apply(g);

      mManager->setGhostsAreStale();
    }

    inline auto getLocalNDHostView() const
    {
      return mManager->getNDHostSubView(mManager->configMemorySizes(), mManager->configLocalSlicing());
    }
    inline auto getFullNDHostView() const { return mManager->getNDHostView(mManager->configMemorySizes()); }
    inline auto getRawHostView() const { return mManager->getRawHostView(); }

    template <typename R> void operator=(R &&g) { this->assign(std::forward<R>(g)); }

    template <typename R> void operator+=(R &&g) { this->operator=(*this + g); }

    void operator=(const ConfigView<T, NDim> &other) { this->assign(other); }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION T &getSet(IDX &&...idx) const
    {
      return mView(idx...);
    }

    inline void confirmSpace(const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      switch (spaceType) {
      case SpaceStateType::Fourier:
        if (!mDisableFFTBlocking)
          throw FieldViewConfigWrongSpaceConfirmation(
              "FieldViewConfig explicitly only can be used in configuration space. Do not transform to Fourier space "
              "in place. Or, if you know what you are doing and you are not doing multiple in-place FFT's on your "
              "integration data, you can call Field<T>::setDisableFFTBlocking() to disable this block, and enable "
              "going from configuration to Fourier space.");
        break;
      case SpaceStateType::Configuration:
      default:
        AbstractField<T, NDim>::confirmSpace(newLayout, spaceType);
        break;
      }
    }

    const LayoutStruct<NDim> &getLayout() const { return mManager->configLayout(); }

    void updateGhosts() const { this->mManager->updateGhosts(); }

    std::string toString() const { return mManager->getName() + "(x)"; }

    /** @brief Disable the blocking of going from configuration space to fourier space.
     *  Used by PowerSpectrumBuilder, which uses newly allocate memory which is filled in configuration space,
     *  and then FFT'ed to fourier space.
     */
    void setDisableFFTBlocking() { mDisableFFTBlocking = true; }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return mView(idx...);
    }

    template <typename R> void onBeforeAssignment(R &&g)
    {
      /* likewise, make sure we are in configuration space (here the FFT may be fired!). */
      mManager->confirmConfigSpace();

      ConfirmSpace::apply(g, mManager->configLayout(), SpaceStateType::Configuration);

      GhostsHunter::apply(g);
      mManager->flagHostMirrorOutdated();
    }

    std::string to_string() const { return mManager->getName() + "(x)"; }

  private:
    // mView is the ONLY member the device kernel reads: eval() above is `return
    // mView(idx...)` and nothing else. Everything else here is host-side bookkeeping that
    // every expression node nevertheless deep-copies, because nodes store their operands by
    // value with the deduced type. mRawView (written in the ctor, never read) and mHostView
    // (declared, never even written) used to sit here and cost 80 B of every Field for
    // nothing -- getRawHostView()/getFullNDHostView() go through mManager->, not through
    // them. mLayout (288 B), memorySizes (24 B) and localSlicing (48 B) followed them into
    // MemoryManager::configLayout() and friends: one copy per field instead of one per copy
    // of every field.
    device::memory::NDViewUnmanaged<T, NDim> mView;

    // Stays here on purpose, and it is the one member that must. It is genuinely per-object
    // mutable state: measuringtools/toolwithownmemory.h sets it on a Field *copy* that shares
    // a manager with a persistent Field, so today the persistent one keeps it false while the
    // copy has it true. In the shared manager it would latch for both and silently weaken the
    // guard in confirmSpace() above. It is 8 bytes.
    bool mDisableFFTBlocking;
  };
} // namespace TempLat

#endif
