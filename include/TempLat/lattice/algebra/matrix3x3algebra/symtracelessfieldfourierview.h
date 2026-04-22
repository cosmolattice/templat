#ifndef TEMPLAT_LATTICE_ALGEBRA_MATRIX3X3ALGEBRA_SYMTRACELESSFIELDFOURIERVIEW_H
#define TEMPLAT_LATTICE_ALGEBRA_MATRIX3X3ALGEBRA_SYMTRACELESSFIELDFOURIERVIEW_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/symtracelessget.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/parallel/device.h"
#include "TempLat/lattice/field/views/fieldviewfourier.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"

namespace TempLat
{

  /** @brief A class which holds symmetric-traceless field in fourier space.
   *
   *
   * Unit test: ctest -R test-complexfieldfourierview
   **/
  template <typename T, size_t NDim> class SymTracelessFieldFourierView
  {
  public:
    // Put public methods here. These should change very little over time.
    SymTracelessFieldFourierView(FourierView<T, NDim> c0, FourierView<T, NDim> c1, FourierView<T, NDim> c2,
                                 FourierView<T, NDim> c3, FourierView<T, NDim> c4)
        : m0(c0), m1(c1), m2(c2), m3(c3), m4(c4), mToolBox(m0.getToolBox()),
          mLayout(mToolBox->mLayouts.getFourierSpaceLayout())
    {
    }

    std::string toString() const
    {
      return "(" + m0.toString() + ", " + m1.toString() + ", " + m2.toString() + ", " + m3.toString() + ", " +
             m4.toString() + ")";
    }

    auto &SymTracelessGet(Tag<0> t) { return m0; }
    const auto &SymTracelessGet(Tag<0> t) const { return m0; }
    auto &operator()(Tag<0> t) { return m0; }
    const auto &operator()(Tag<0> t) const { return m0; }

    auto &SymTracelessGet(Tag<1> t) { return m1; }
    const auto &SymTracelessGet(Tag<1> t) const { return m1; }
    auto &operator()(Tag<1> t) { return m1; }
    const auto &operator()(Tag<1> t) const { return m1; }

    auto &SymTracelessGet(Tag<2> t) { return m2; }
    const auto &SymTracelessGet(Tag<2> t) const { return m2; }
    auto &operator()(Tag<2> t) { return m2; }
    const auto &operator()(Tag<2> t) const { return m2; }

    auto &SymTracelessGet(Tag<3> t) { return m3; }
    const auto &SymTracelessGet(Tag<3> t) const { return m3; }
    auto &operator()(Tag<3> t) { return m3; }
    const auto &operator()(Tag<3> t) const { return m3; }

    auto &SymTracelessGet(Tag<4> t) { return m4; }
    const auto &SymTracelessGet(Tag<4> t) const { return m4; }
    auto &operator()(Tag<4> t) { return m4; }
    const auto &operator()(Tag<4> t) const { return m4; }

    auto &SymTracelessGet(Tag<1> t1, Tag<1> t2) { return m0; }
    const auto &SymTracelessGet(Tag<1> t1, Tag<1> t2) const { return m0; }
    auto &operator()(Tag<1> t1, Tag<1> t2) { return m0; }
    const auto &operator()(Tag<1> t1, Tag<1> t2) const { return m0; }

    auto &SymTracelessGet(Tag<1> t1, Tag<2> t2) { return m1; }
    const auto &SymTracelessGet(Tag<1> t1, Tag<2> t2) const { return m1; }
    auto &operator()(Tag<1> t1, Tag<2> t2) { return m1; }
    const auto &operator()(Tag<1> t1, Tag<2> t2) const { return m1; }

    auto &SymTracelessGet(Tag<1> t1, Tag<3> t2) { return m2; }
    const auto &SymTracelessGet(Tag<1> t1, Tag<3> t2) const { return m2; }
    auto &operator()(Tag<1> t1, Tag<3> t2) { return m2; }
    const auto &operator()(Tag<1> t1, Tag<3> t2) const { return m2; }

    auto &SymTracelessGet(Tag<2> t1, Tag<1> t2) { return m1; }
    const auto &SymTracelessGet(Tag<2> t1, Tag<1> t2) const { return m1; }
    auto &operator()(Tag<2> t1, Tag<1> t2) { return m1; }
    const auto &operator()(Tag<2> t1, Tag<1> t2) const { return m1; }

    auto &SymTracelessGet(Tag<2> t1, Tag<2> t2) { return m3; }
    const auto &SymTracelessGet(Tag<2> t1, Tag<2> t2) const { return m3; }
    auto &operator()(Tag<2> t1, Tag<2> t2) { return m3; }
    const auto &operator()(Tag<2> t1, Tag<2> t2) const { return m3; }

    auto &SymTracelessGet(Tag<2> t1, Tag<3> t2) { return m4; }
    const auto &SymTracelessGet(Tag<2> t1, Tag<3> t2) const { return m4; }
    auto &operator()(Tag<2> t1, Tag<3> t2) { return m4; }
    const auto &operator()(Tag<2> t1, Tag<3> t2) const { return m4; }

    auto &SymTracelessGet(Tag<3> t1, Tag<1> t2) { return m2; }
    const auto &SymTracelessGet(Tag<3> t1, Tag<1> t2) const { return m2; }
    auto &operator()(Tag<3> t1, Tag<1> t2) { return m2; }
    const auto &operator()(Tag<3> t1, Tag<1> t2) const { return m2; }

    auto &SymTracelessGet(Tag<3> t1, Tag<2> t2) { return m4; }
    const auto &SymTracelessGet(Tag<3> t1, Tag<2> t2) const { return m4; }
    auto &operator()(Tag<3> t1, Tag<2> t2) { return m4; }
    const auto &operator()(Tag<3> t1, Tag<2> t2) const { return m4; }

    auto SymTracelessGet(Tag<3> t1, Tag<3> t2)
    {
      return -m0 - m3;
    } // TODO: Jorge: I have been forced to remove the & to do assignements. I do not like this, so we should discuss
      // about this.
    const auto SymTracelessGet(Tag<3> t1, Tag<3> t2) const { return -m0 - m3; }
    auto operator()(Tag<3> t1, Tag<3> t2) { return -m0 - m3; }
    const auto operator()(Tag<3> t1, Tag<3> t2) const { return -m0 - m3; }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<decltype(m0.eval(idx...)), 5> result;
      result[0] = m0.eval(idx...);
      result[1] = m1.eval(idx...);
      result[2] = m2.eval(idx...);
      result[3] = m3.eval(idx...);
      result[4] = m4.eval(idx...);
      return result;
    }

    template <typename R> void operator=(R &&g)
    {
      const auto &g0 = SymTracelessGetter::get(g, 0_c);
      const auto &g1 = SymTracelessGetter::get(g, 1_c);
      const auto &g2 = SymTracelessGetter::get(g, 2_c);
      const auto &g3 = SymTracelessGetter::get(g, 3_c);
      const auto &g4 = SymTracelessGetter::get(g, 4_c);

      m0.onBeforeAssignment(g0);
      m1.onBeforeAssignment(g1);
      m2.onBeforeAssignment(g2);
      m3.onBeforeAssignment(g3);
      m4.onBeforeAssignment(g4);

      PreGet::apply(g);

      const auto view0 = m0.getView();
      const auto view1 = m1.getView();
      const auto view2 = m2.getView();
      const auto view3 = m3.getView();
      const auto view4 = m4.getView();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply(
            [&](auto &&...args) {
              auto result = DoEval::eval(g, args...);
              view0(args...) = result[0];
              view1(args...) = result[1];
              view2(args...) = result[2];
              view3(args...) = result[3];
              view4(args...) = result[4];
            },
            idx);
      };
      device::iteration::foreach ("SymTracelessFourierViewAssign", mLayout, functor);

      PostGet::apply(g);
    }

    auto getDx() const { return m0.getDx(); }
    auto getKIR() const { return m0.getKIR(); }

    using Getter = SymTracelessGetter;
    static constexpr size_t SHIFTIND = 0;
    static constexpr size_t size = 5;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    FourierView<T, NDim> m0;
    FourierView<T, NDim> m1;
    FourierView<T, NDim> m2;
    FourierView<T, NDim> m3;
    FourierView<T, NDim> m4;

    device::memory::host_ptr<MemoryToolBox<NDim>> mToolBox;

    LayoutStruct<NDim> mLayout;
  };
} // namespace TempLat

#endif
