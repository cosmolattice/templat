#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSFIELD_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSFIELD_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/parallel/device.h"
#include "TempLat/lattice/field/assignablefieldcollection.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/symtracelessget.h"
#include "TempLat/lattice/algebra/helpers/getdx.h"
#include "TempLat/lattice/algebra/helpers/getkir.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelesswrapper.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfieldfourierview.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include <memory>

namespace TempLat
{
  /** @brief A class which implements symmetric-traceless matrix field used for GWs.
   *
   * Unit test: ctest -R test-complexfield
   **/
  template <typename T, size_t _NDim = 0> class SymTracelessField
  {
  public:
    // Put public methods here. These should change very little over time.

    static_assert(_NDim != 0, "NDim template parameter is required. Use e.g. SymTracelessField<double, 3>.");

    static constexpr size_t NDim = _NDim;

    SymTracelessField(Field<T, NDim> f0, Field<T, NDim> f1, Field<T, NDim> f2, Field<T, NDim> f3, Field<T, NDim> f4):
    m0(f0),
    m1(f1),
    m2(f2),
    m3(f3),
    m4(f4),
    mName("symtraceless(" + f0.getName() + ", " + f1.getName() + ", " + f2.getName() + ", " + f3.getName() + ", " + f4.getName() + ")"),
    mToolBox(m1.getToolBox() == nullptr ? m1.getToolBox() : m0.getToolBox()),
    mLayout(mToolBox->mLayouts.getConfigSpaceLayout())
    {
    }

    SymTracelessField(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
                 LatticeParameters<T> pLatPar = LatticeParameters<T>()):
    m0("c11_" + name, toolBox, pLatPar),
    m1("c12_" + name, toolBox, pLatPar),
    m2("c13_" + name, toolBox, pLatPar),
    m3("c22_" + name, toolBox, pLatPar),
    m4("c23_" + name, toolBox, pLatPar),
    mName(name), mToolBox(toolBox),
    mLayout(mToolBox->mLayouts.getConfigSpaceLayout())
    {
    }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<0> t) { return m0; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<0> t) const { return m0; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<0> t) { return m0; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<0> t) const { return m0; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<1> t) { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<1> t) const { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<1> t) { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<1> t) const { return m1; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<2> t) { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<2> t) const { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<2> t) { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<2> t) const { return m2; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<3> t) { return m3; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<3> t) const { return m3; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<3> t) { return m3; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<3> t) const { return m3; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<4> t) { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<4> t) const { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<4> t) { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<4> t) const { return m4; }


    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<1> t1, Tag<1> t2) { return m0; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<1> t1, Tag<1> t2) const { return m0; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<1> t1, Tag<1> t2) { return m0; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<1> t1, Tag<1> t2) const { return m0; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<1> t1, Tag<2> t2) { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<1> t1, Tag<2> t2) const { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<1> t1, Tag<2> t2) { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<1> t1, Tag<2> t2) const { return m1; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<1> t1, Tag<3> t2) { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<1> t1, Tag<3> t2) const { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<1> t1, Tag<3> t2) { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<1> t1, Tag<3> t2) const { return m2; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<2> t1, Tag<1> t2) { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<2> t1, Tag<1> t2) const { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<2> t1, Tag<1> t2) { return m1; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<2> t1, Tag<1> t2) const { return m1; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<2> t1, Tag<2> t2) { return m3; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<2> t1, Tag<2> t2) const { return m3; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<2> t1, Tag<2> t2) { return m3; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<2> t1, Tag<2> t2) const { return m3; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<2> t1, Tag<3> t2) { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<2> t1, Tag<3> t2) const { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<2> t1, Tag<3> t2) { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<2> t1, Tag<3> t2) const { return m4; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<3> t1, Tag<1> t2) { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<3> t1, Tag<1> t2) const { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<3> t1, Tag<1> t2) { return m2; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<3> t1, Tag<1> t2) const { return m2; }

    DEVICE_FORCEINLINE_FUNCTION
    auto &SymTracelessGet(Tag<3> t1, Tag<2> t2) { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &SymTracelessGet(Tag<3> t1, Tag<2> t2) const { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    auto &operator()(Tag<3> t1, Tag<2> t2) { return m4; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto &operator()(Tag<3> t1, Tag<2> t2) const { return m4; }

    DEVICE_FORCEINLINE_FUNCTION
    auto SymTracelessGet(Tag<3> t1, Tag<3> t2) { return - m0 - m3; } //TODO: Jorge: I have been forced to remove the & to do assignements. I do not like this, so we should discuss about this.
    DEVICE_FORCEINLINE_FUNCTION
    const auto SymTracelessGet(Tag<3> t1, Tag<3> t2) const { return - m0 - m3; }
    DEVICE_FORCEINLINE_FUNCTION
    auto operator()(Tag<3> t1, Tag<3> t2) { return - m0 - m3; }
    DEVICE_FORCEINLINE_FUNCTION
    const auto operator()(Tag<3> t1, Tag<3> t2) const { return - m0 - m3; }


    template <int N> auto &operator()(Tag<N> t) { return SymTracelessGet(t); }
    template <int N> const auto &operator()(Tag<N> t) const { return SymTracelessGet(t); }

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<T, 5> result;
      result[0] = m0.eval(idx...);
      result[1] = m1.eval(idx...);
      result[2] = m2.eval(idx...);
      result[3] = m3.eval(idx...);
      result[4] = m4.eval(idx...);
      return result;
    }

    SymTracelessFieldFourierView<T, NDim> inFourierSpace() { return {m0.inFourierSpace(),
                                                                           m1.inFourierSpace(),
                                                                           m2.inFourierSpace(),
                                                                           m3.inFourierSpace(),
                                                                           m4.inFourierSpace()}; }

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
      device::iteration::foreach ("SymTracelessConfigViewAssign", mLayout, functor);

      PostGet::apply(g);

      m0.setGhostsAreStale();
      m1.setGhostsAreStale();
      m2.setGhostsAreStale();
      m3.setGhostsAreStale();
      m4.setGhostsAreStale();
    }

    template <typename R> void operator+=(R &&r) { (*this) = (*this) + r; }

    std::string toString() const { return *mName; }

    DEVICE_FORCEINLINE_FUNCTION
    auto getDx() const { return m0.getDx(); }

    DEVICE_FORCEINLINE_FUNCTION
    auto getKIR() const { return m0.getKIR(); }

    void updateGhosts()
    {
      m0.updateGhosts();
      m1.updateGhosts();
      m2.updateGhosts();
      m3.updateGhosts();
      m4.updateGhosts();
    }

    void setGhostsAreStale()
    {
      m0.setGhostsAreStale();
      m1.setGhostsAreStale();
      m2.setGhostsAreStale();
      m3.setGhostsAreStale();
      m4.setGhostsAreStale();
    }

    using Getter = SymTracelessGetter;
    static constexpr size_t SHIFTIND = 0;
    static constexpr size_t size = 5;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    Field<T, NDim> m0;
    Field<T, NDim> m1;
    Field<T, NDim> m2;
    Field<T, NDim> m3;
    Field<T, NDim> m4;

    device::memory::host_string mName;

    device::memory::host_ptr<MemoryToolBox<NDim>> mToolBox;

    LayoutStruct<NDim> mLayout;
  };

} // namespace TempLat

#endif
