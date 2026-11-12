#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELD_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELD_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#include "TempLat/parallel/device.h"
#include "TempLat/lattice/field/assignablefieldcollection.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/complexfieldget.h"
#include "TempLat/lattice/algebra/helpers/getdx.h"
#include "TempLat/lattice/algebra/helpers/getkir.h"
#include "TempLat/lattice/algebra/complexalgebra/complexwrapper.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldfourierview.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/confirmghosts.h"
#include <memory>

namespace TempLat
{
  /** @brief A class which implements complex-valued fields.
   *
   * Unit test: ctest -R test-complexfield
   *
   * @vocab-summary A complex scalar field, stored as two real components.
   * @vocab-signature ComplexField<T, NDim> psi("psi", toolBox);
   **/
  template <typename T, size_t _NDim = 0> class ComplexField
  {
  public:
    // Put public methods here. These should change very little over time.
    static_assert(_NDim != 0, "NDim template parameter is required. Use e.g. ComplexField<double, 3>.");

    static constexpr size_t NDim = _NDim;

    ComplexField(Field<T, NDim> f1, Field<T, NDim> f2)
        : mR(f1), mI(f2),
          mName("complex(" + f1.getMemoryManager()->getName() + ", " + f2.getMemoryManager()->getName() + ")"),
          mToolBox(mR.getToolBox() == nullptr ? mI.getToolBox() : mR.getToolBox()),
          mLayout(mToolBox->mLayouts.getConfigSpaceLayout())
    {
    }

    ComplexField(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
                 LatticeParameters<T> pLatPar = LatticeParameters<T>())
        : mR("Re_" + name, toolBox, pLatPar), mI("Im_" + name, toolBox, pLatPar), mName(name), mToolBox(toolBox),
          mLayout(mToolBox->mLayouts.getConfigSpaceLayout())
    {
    }

    auto &ComplexFieldGet(Tag<0> t) { return mR; }
    const auto &ComplexFieldGet(Tag<0> t) const { return mR; }
    auto &operator()(Tag<0> t) { return mR; }
    const auto &operator()(Tag<0> t) const { return mR; }

    auto &ComplexFieldGet(Tag<1> t) { return mI; }
    const auto &ComplexFieldGet(Tag<1> t) const { return mI; }
    auto &operator()(Tag<1> t) { return mI; }
    const auto &operator()(Tag<1> t) const { return mI; }

    template <int N> auto &operator()(Tag<N> t)
    {
      static_assert(N >= 0 && N <= 1, "Operator(): N must be 0 or 1 for ComplexField");
      return ComplexFieldGet(t);
    }
    template <int N> const auto &operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 1, "Operator(): N must be 0 or 1 for ComplexField");
      return ComplexFieldGet(t);
    }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<T, 2> result;
      result[0] = mR.eval(idx...);
      result[1] = mI.eval(idx...);
      return result;
    }

    ComplexFieldFourierView<T, NDim> inFourierSpace() { return {mR.inFourierSpace(), mI.inFourierSpace()}; }

    /** @brief The fields an assignment writes; they receive components 0 and 1 of the value. Lets collections assign
     *  all members in one kernel, see FieldCollection::operator=. */
    auto assignedFields() const { return device::array<Field<T, NDim>, 2>{mR, mI}; }
    static constexpr size_t firstAssignedComponent = 0;

    /** @brief Prepares the right-hand side g for the single assignment kernel writing all components, see
     *  prepareAssignment. */
    template <typename R> GhostReach onBeforeAssignment(R &&g)
    {
      return prepareAssignment<T, NDim>(g, componentManagers(), mLayout);
    }

    template <typename R> void operator=(R &&g)
    {
      const GhostReach targetReach = onBeforeAssignment(g);
      assignViews("ComplexConfigViewAssign", device::array<decltype(mR.getView()), 2>{mR.getView(), mI.getView()}, g,
                  mLayout, targetReach, mR.toolBox()->mSweepBuffer);

      mR.setGhostsAreStale();
      mI.setGhostsAreStale();
    }

    template <typename R> void operator+=(R &&r) { (*this) = (*this) + r; }

    std::string toString() const { return *mName; }

    /** @brief Ghost walk entry batched over the components, see ConfirmGhosts. Without a reach, updates the
     *  ghosts of all components if they are stale. */
    device::Idx confirmGhostsUpToDate(const GhostReach &reach) const
    {
      return MemoryManager<T, NDim>::confirmGhostsUpToDateBatch(componentManagers(), reach);
    }
    device::Idx confirmGhostsUpToDate() const
    {
      return MemoryManager<T, NDim>::confirmGhostsUpToDateBatch(componentManagers());
    }
    std::array<MemoryManager<T, NDim> *, 2> componentManagers() const
    {
      return {mR.getMemoryManager().get(), mI.getMemoryManager().get()};
    }
    void confirmSpace(const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      ConfirmSpace::apply(mR, newLayout, spaceType);
      ConfirmSpace::apply(mI, newLayout, spaceType);
    }

    auto getDx() const { return mR.getDx(); }
    auto getKIR() const { return mR.getKIR(); }

    void updateGhosts() { MemoryManager<T, NDim>::updateGhostsBatch(componentManagers()); }

    void setGhostsAreStale()
    {
      mR.setGhostsAreStale();
      mI.setGhostsAreStale();
    }

    using Getter = ComplexFieldGetter;
    static constexpr size_t SHIFTIND = 0;
    static constexpr size_t size = 2;

    device::memory::host_ptr<MemoryToolBox<NDim>> getToolBox() const { return mToolBox; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    Field<T, NDim> mR;
    Field<T, NDim> mI;

    device::memory::host_string mName;

    device::memory::host_ptr<MemoryToolBox<NDim>> mToolBox;

    LayoutStruct<NDim> mLayout;
  };

  /**
   * @vocab-summary Bundles two existing real Fields into one complex field, sharing their memory rather than
   * copying it.
   * @vocab-signature CField(re, im)
   **/
  template <typename T, size_t NDim> auto CField(Field<T, NDim> f1, Field<T, NDim> f2)
  {
    return ComplexField<T, NDim>(f1, f2);
  }
} // namespace TempLat

#endif
