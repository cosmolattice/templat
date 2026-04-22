#ifndef TEMPLAT_LATTICE_ALGEBRA_COMPLEXALGEBRA_COMPLEXFIELDFOURIERVIEW_H
#define TEMPLAT_LATTICE_ALGEBRA_COMPLEXALGEBRA_COMPLEXFIELDFOURIERVIEW_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/complexalgebra/helpers/complexfieldget.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/parallel/device.h"
#include "TempLat/lattice/field/views/fieldviewfourier.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"

namespace TempLat
{

  /** @brief A class which holds complex field in fourier space.
   *
   *
   * Unit test: ctest -R test-complexfieldfourierview
   **/
  template <typename T, size_t NDim> class ComplexFieldFourierView
  {
  public:
    // Put public methods here. These should change very little over time.
    ComplexFieldFourierView(FourierView<T, NDim> r, FourierView<T, NDim> i)
        : mR(r), mI(i), mToolBox(mR.getToolBox() == nullptr ? mI.getToolBox() : mR.getToolBox()),
          mLayout(mToolBox->mLayouts.getFourierSpaceLayout())
    {
    }

    std::string toString() const { return "(" + mR.toString() + ", " + mI.toString() + ")"; }

    auto &ComplexFieldGet(Tag<0> t) { return mR; }
    const auto &ComplexFieldGet(Tag<0> t) const { return mR; }
    auto &operator()(Tag<0> t) { return mR; }
    const auto &operator()(Tag<0> t) const { return mR; }

    auto &ComplexFieldGet(Tag<1> t) { return mI; }
    const auto &ComplexFieldGet(Tag<1> t) const { return mI; }
    auto &operator()(Tag<1> t) { return mI; }
    const auto &operator()(Tag<1> t) const { return mI; }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<decltype(mR.eval(idx...)), 2> result;
      result[0] = mR.eval(idx...);
      result[1] = mI.eval(idx...);
      return result;
    }

    template <typename R> void operator=(R &&g)
    {
      const auto &gR = g.ComplexFieldGet(0_c);
      const auto &gI = g.ComplexFieldGet(1_c);

      mR.onBeforeAssignment(gR);
      mI.onBeforeAssignment(gI);

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
      device::iteration::foreach ("ComplexFourierViewAssign", mLayout, functor);

      PostGet::apply(g);
    }

    auto getDx() const { return mR.getDx(); }
    auto getKIR() const { return mR.getKIR(); }

    using Getter = ComplexFieldGetter;
    static constexpr size_t SHIFTIND = 0;
    static constexpr size_t size = 2;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    FourierView<T, NDim> mR;
    FourierView<T, NDim> mI;

    device::memory::host_ptr<MemoryToolBox<NDim>> mToolBox;

    LayoutStruct<NDim> mLayout;
  };
} // namespace TempLat

#endif
