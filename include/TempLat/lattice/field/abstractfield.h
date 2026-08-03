#ifndef TEMPLAT_LATTICE_FIELD_ABSTRACTFIELD_H
#define TEMPLAT_LATTICE_FIELD_ABSTRACTFIELD_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/spacestateinterface.h"
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/latticeparameters.h"
#include "TempLat/lattice/memory/memorylayoutstate.h"
#include "TempLat/lattice/memory/memorymanager.h"
#include "TempLat/lattice/memory/memorytoolbox.h"

#include "TempLat/lattice/algebra/helpers/getderiv.h"

#ifdef TEMPLAT_FIELD_COPY_COUNTER
#include "TempLat/parallel/device.h"
#include <atomic>
#endif

namespace TempLat
{

  MakeException(FieldValueGetterException);

  enum CANONICALTYPE { AMPLITUDE, MOMENTUM };

  /** @brief A base class for a one-component field, providing common functionality and holding relevant pointers to
   * memory tools.
   *
   **/
  template <typename T, size_t _NDim> class AbstractField
  {
  public:
    // Put public methods here. These should change very little over time.

    static constexpr size_t NDim = _NDim;

    AbstractField(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, LatticeParameters<T> pLatPar)
        : mManager(toolBox, name), latPar(pLatPar)
    {
    }

#ifdef TEMPLAT_FIELD_COPY_COUNTER
    /** @brief Measurement instrument -- OFF unless -DTEMPLAT_FIELD_COPY_COUNTER.
     *
     * Expression-template nodes store their operands by value, so building one expression
     * deep-copies every leaf once per level of tree depth. That multiplier had only ever
     * been inferred from timing; this counts it.
     *
     * AbstractField is a subobject of BOTH views of a Field (ConfigView's base and the
     * embedded FourierView's base), so ONE Field copy raises this by TWO.
     *
     * The atomic RMW costs real time, so a binary built with this macro must not be used
     * for timings -- build a separate, clean one for those. With the macro undefined this
     * block disappears entirely and the build is unchanged.
     */
    static inline std::atomic<long> sCopyCount{0};

    // No mToolBox here: this branch reaches the toolbox through the manager, so the
    // member is gone. Keep the two init lists in step when this merges either way.
    DEVICE_INLINE_FUNCTION AbstractField(const AbstractField &o) : mManager(o.mManager), latPar(o.latPar)
    {
      // The counter is host state, but the copy constructor must stay callable from
      // __host__ __device__ context: nvcc instantiates it wherever the implicit one was
      // reachable before.
      KOKKOS_IF_ON_HOST((sCopyCount.fetch_add(1, std::memory_order_relaxed);))
    }

    // Declaring a copy constructor suppresses the implicit move constructor, which would
    // silently route every move through the (now counting, and more expensive) copy. Put it
    // back explicitly so the instrumented build moves exactly where the clean build does.
    AbstractField(AbstractField &&) = default;
#endif

    inline void confirmSpace(const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      switch (spaceType) {
      case SpaceStateType::Fourier:
        mManager->confirmFourierSpace();
        break;
      case SpaceStateType::Configuration:
      default:
        mManager->confirmConfigSpace();
        break;
      }
    }

    device::memory::host_ptr<MemoryToolBox<NDim>> getToolBox() const { return mManager->getToolBox(); }

    device::Idx confirmGhostsUpToDate() const { return this->mManager->confirmGhostsUpToDate(); }

    // Mostly for testing purpose

    /** @brief Check the current state.
     *
     * @return true if in configuration space
     */
    bool isConfigSpace() const { return mManager->isConfigSpace(); }

    /** @brief Check the current state.
     *
     * @return true if in fourier space
     */
    bool isFourierSpace() const { return mManager->isFourierSpace(); }

    void setGhostsAreStale() const { mManager->setGhostsAreStale(); }
    bool areGhostsStale() const { return mManager->areGhostsStale(); }

    device::memory::host_ptr<MemoryManager<T, NDim>> getMemoryManager() const { return mManager; }

    const BCSpec<NDim> &getBCSpec() const { return mManager->getBCSpec(); }

    void setBCSpec(BCSpec<NDim> bcSpec)
    {
      mManager->setBCSpec(bcSpec);
      mManager->setGhostsAreStale();
    }

    auto getDx() const { return latPar.getDx(); }
    auto getKIR() const { return latPar.getKIR(); }
    auto getLatParams() const { return latPar; }

  protected:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    /** @brief The toolbox, without copying a reference count.
     *
     * Replaces the former mToolBox member: the manager already holds the toolbox
     * (memorymanager.h), so a field-side handle was a second reference count paid on
     * every copy of a field -- and fields are copied by value into every node of every
     * expression template. Use this in field internals; getToolBox() above still hands
     * out an owning copy for external callers.
     */
    const device::memory::host_ptr<MemoryToolBox<NDim>> &toolBox() const { return mManager->getToolBox(); }

    device::memory::host_ptr<MemoryManager<T, NDim>> mManager;

    const LatticeParameters<T> latPar; // Information about the lattice (dx, kir...)
                                       // Conceptually not amazing but really useful.
  };
} // namespace TempLat

#endif
