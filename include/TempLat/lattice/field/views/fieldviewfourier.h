#ifndef TEMPLAT_LATTICE_FIELD_VIEWS_FIELDVIEWFOURIER_H
#define TEMPLAT_LATTICE_FIELD_VIEWS_FIELDVIEWFOURIER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/lattice/field/abstractfield.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/parallel/device_iteration.h"

namespace TempLat
{
  MakeException(FieldViewFourierWrongSpaceConfirmation);

  /** @brief A view on the field which, when interacted with, assures every time again that things are in
   *   *fourier* space, and possibly the ghost cells are updated when needed.
   *
   *   No public constructor: only the friend class Field can instantiate this view on its own memory.
   *
   * Unit test: ctest -R test-fieldviewconfig
   **/
  template <typename T, size_t _NDim> class FourierView : public AbstractField<T, _NDim>
  {
  public:
    static constexpr size_t NDim = _NDim;

    using AbstractField<T, NDim>::mManager;
    using AbstractField<T, NDim>::mToolBox;

    template <typename R> void operator=(R &&g) { this->assign(std::forward<R>(g)); }

    template <typename R> void assign(R &&g)
    {
      const auto layout = mToolBox->mLayouts.getFourierSpaceLayout();
      onBeforeAssignment(g);

      PreGet::apply(g);

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply([&](auto &&...args) { mView(args...) = DoEval::eval(g, args...); }, idx);
      };
      device::iteration::foreach ("FourierViewAssign", layout, functor);

      PostGet::apply(g);
    }

    template <typename R>
      requires requires(R r) {
        r.ComplexFieldGet(0_c);
        r.ComplexFieldGet(1_c);
      }
    void assign(R &&g)
    {
      const auto layout = mToolBox->mLayouts.getFourierSpaceLayout();
      onBeforeAssignment(g);

      PreGet::apply(g);

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply(
            [&](auto &&...args) {
              const auto result = DoEval::eval(g, args...);
              mView(args...) = complex<T>(result[0], result[1]);
            },
            idx);
      };
      device::iteration::foreach ("FourierViewAssign", layout, functor);

      PostGet::apply(g);
    }

    inline auto getLocalNDHostView() const
    {
      // As we have no ghosts in Fourier space, we can return the full view.
      return mManager->template getNDHostView<complex<T>>(memorySizes);
    }
    inline auto getFullNDHostView() const { return mManager->template getNDHostView<complex<T>>(memorySizes); }
    inline auto getRawHostView() const { return mManager->template getRawHostView<complex<T>>(); }

    template <typename R> void onBeforeAssignment(R &&g)
    {
      /* likewise, make sure we are in configuration space (here the FFT may be fired!). */
      mManager->confirmFourierSpace();
      ConfirmSpace::apply(g, mToolBox->mLayouts.getFourierSpaceLayout(), SpaceStateType::Fourier);
      GhostsHunter::apply(g);
      mManager->flagHostMirrorOutdated();
    }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION complex<T> &getSet(const IDX &...idx) const
    {
      return mView(idx...);
    }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return mView(idx...);
    }

    inline void confirmSpace(const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      switch (spaceType) {
      case SpaceStateType::Configuration:
        throw FieldViewFourierWrongSpaceConfirmation("FieldViewFourier explicitly only can be used in Fourier space. "
                                                     "Do not transform to configuration space in place.");
        break;
      case SpaceStateType::Fourier:
      default:
        AbstractField<T, NDim>::confirmSpace(newLayout, spaceType);
        break;
      }
    }

    auto getView() const { return mView; }

    std::string toString() const { return mManager->getName() + "(k)"; }

    const auto &getLayout() const { return mToolBox->mLayouts.getFourierSpaceLayout(); }

    // MPI aware setting of value. Use exceptionally (remove zero mode for example)

    template <typename... Args> void setZeroMode(const complex<T> &toSet)
    {
      // This is dimension-aware.
      device::IdxArray<NDim> global_coord{{}};
      device::IdxArray<NDim> mem_pos{{}};

      const auto &layout = mToolBox->mLayouts.getFourierSpaceLayout();
      const bool owned = device::apply(
          [&](const auto &...idx) { return layout.putMemoryIndexFromSpatialLocationInto(mem_pos, idx...); },
          global_coord);

      // do this only if this process owns the zero mode!
      if (owned) device::memory::setAtOnePoint(*this, mem_pos, toSet);
    }

    std::string to_string() const { return mManager->getName() + "(k)"; }

    template <typename S, size_t __NDim> friend class Field;

  private:
    FourierView(const AbstractField<T, NDim> &f) : AbstractField<T, NDim>(f)
    {
      if (mToolBox == nullptr) return;
      auto layout = mToolBox->mLayouts.getFourierSpaceLayout();
      memorySizes = layout.getSizesInMemory();
      mView = mManager->template getNDView<complex<T>>(memorySizes);
    }

    // Same story as ConfigView, and it costs twice over: every Field embeds one of these,
    // AND assign() above captures [=, *this], so the dead members were also in the kernel
    // argument capture of every Fourier-space assignment. mRawView (written, never read),
    // mHostView and localSlicing (declared, never even written) were 128 B of every Field.
    // memorySizes IS live, in getLocalNDHostView/getFullNDHostView above.
    device::memory::NDViewUnmanaged<complex<T>, NDim> mView;

    device::IdxArray<NDim> memorySizes;
  };
} // namespace TempLat

#endif
