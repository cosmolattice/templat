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

    /** @brief Confirms space and ghosts for the right-hand side g. Returns true if g reads one of this field's
     *  components at a neighbouring site, i.e. if assigning g in place would read already overwritten values. */
    template <typename R> bool onBeforeAssignment(R &&g)
    {
      // Walk the whole expression once per component: g reads mI at shifted sites as much as mR.
      const bool readR = mR.onBeforeAssignment(g);
      const bool readI = mI.onBeforeAssignment(g);
      return readR || readI;
    }

    template <typename R> void operator=(R &&g)
    {
      onBeforeAssignment(g);

      PreGet::apply(g);

      const auto viewR = mR.getView();
      const auto viewI = mI.getView();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply(
            [&](auto &&...args) {
              auto result = DoEval::eval(g, args...);
              viewR(args...) = result[0];
              viewI(args...) = result[1];
            },
            idx);
      };
      device::iteration::foreach ("ComplexConfigViewAssign", mLayout, functor);

      PostGet::apply(g);

      mR.setGhostsAreStale();
      mI.setGhostsAreStale();
    }

    template <typename R> void operator+=(R &&r) { (*this) = (*this) + r; }

    std::string toString() const { return *mName; }

    // Space/ghost confirmation of the components, so that this leaf can be walked as part of a whole expression.
    // A bare field needs no ghosts by itself; a shift above it calls confirmGhostsUpToDate.
    void doWeNeedGhosts() const {}
    device::Idx confirmGhostsUpToDate() const
    {
      MemoryManager<T, NDim> *mgrs[] = {mR.getMemoryManager().get(), mI.getMemoryManager().get()};
      return MemoryManager<T, NDim>::confirmGhostsUpToDateBatch(mgrs);
    }
    void confirmSpace(const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      ConfirmSpace::apply(mR, newLayout, spaceType);
      ConfirmSpace::apply(mI, newLayout, spaceType);
    }

    auto getDx() const { return mR.getDx(); }
    auto getKIR() const { return mR.getKIR(); }

    void updateGhosts()
    {
      MemoryManager<T, NDim> *mgrs[] = {mR.getMemoryManager().get(), mI.getMemoryManager().get()};
      MemoryManager<T, NDim>::updateGhostsBatch(mgrs);
    }

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
