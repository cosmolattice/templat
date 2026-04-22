#ifndef COSMOINTERFACE_SU2ALGEBRA_SU2LIEALGEBRAFIELD_H
#define COSMOINTERFACE_SU2ALGEBRA_SU2LIEALGEBRAFIELD_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2021

#include "TempLat/lattice/algebra/su2algebra/su2field.h"

namespace TempLat
{
  /** @brief A class which implements a SU2 field (algebra). Same than the group, except set 0 component to 0 and has a
   *function which returns its components wrt to sigma/2 instead of sigma, which is more conventional for the algebra
   *(but internally, also expanded as a function of sigma).
   *
   * Unit test: ctest -R test-su2liealgebrafield
   **/
  template <typename T, size_t _NDim = 0> class SU2LieAlgebraField
  {
  public:
    // Put public methods here. These should change very little over time.
    static_assert(_NDim != 0, "NDim template parameter is required. Use e.g. SU2LieAlgebraField<double, 3>.");

    static constexpr size_t NDim = _NDim;

    SU2LieAlgebraField(Field<T, NDim> f1, Field<T, NDim> f2, Field<T, NDim> f3)
        : fs{{f1, f2, f3}}, mName("NoName"), mLayout(fs[0].getToolBox()->mLayouts.getConfigSpaceLayout())
    {
    }

    SU2LieAlgebraField(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
                       LatticeParameters<T> pLatPar)
        : fs{{
              Field<T, NDim>(name + "_1", toolBox, pLatPar), //
              Field<T, NDim>(name + "_2", toolBox, pLatPar), //
              Field<T, NDim>(name + "_3", toolBox, pLatPar)  //
          }},
          mName(name), mLayout(toolBox->mLayouts.getConfigSpaceLayout())
    {
    }

    template <int N> auto SU2Get(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "SU2Get: N must be between 0 and 3 for SU2LieAlgebraField");
      return operator()(t);
    }

    template <int N> auto SU2LieAlgebraGet(Tag<N> t) const { return 2 * SU2Get(t); }

    ZeroType operator()(Tag<0> t) const { return {}; }

    template <int M> auto &operator()(Tag<M> t)
    {
      static_assert(M >= 1 && M <= 3, "Operator(): M must be between 0 and 3 for SU2LieAlgebraField");
      return fs[M - 1];
    }
    template <int M> const auto &operator()(Tag<M> t) const
    {
      static_assert(M >= 1 && M <= 3, "Operator(): M must be between 0 and 3 for SU2LieAlgebraField");
      return fs[M - 1];
    }

    template <typename R> void operator=(R &&r)
    {
      fs[0].onBeforeAssignment(r.SU2Get(1_c));
      fs[1].onBeforeAssignment(r.SU2Get(2_c));
      fs[2].onBeforeAssignment(r.SU2Get(3_c));

      PreGet::apply(r);

      const auto view1 = fs[0].getView();
      const auto view2 = fs[1].getView();
      const auto view3 = fs[2].getView();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply(
            [&](const auto &...args) {
              auto result = DoEval::eval(r, args...);
              view1(args...) = result[1];
              view2(args...) = result[2];
              view3(args...) = result[3];
            },
            idx);
      };
      device::iteration::foreach ("SU2AlgebraConfigViewAssign", mLayout, functor);

      PostGet::apply(r);

      fs[0].setGhostsAreStale();
      fs[1].setGhostsAreStale();
      fs[2].setGhostsAreStale();
    }

    template <typename... IDX>
      requires requires(Field<T, NDim> f, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(f, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<T, 4> result;
      result[0] = T(0);
      result[1] = fs[0].eval(idx...);
      result[2] = fs[1].eval(idx...);
      result[3] = fs[2].eval(idx...);
      return result;
    }

    std::string toString() const { return *mName; }

    auto getDx() const { return GetDx::getDx(fs[0]); }
    auto getKIR() const { return GetKIR::getKIR(fs[0]); }

    inline auto getToolBox() { return GetToolBox::get(fs[0]); }

    inline void updateGhosts()
    {
      fs[0].updateGhosts();
      fs[1].updateGhosts();
      fs[2].updateGhosts();
    }

    using Getter = SU2Getter;
    static constexpr size_t SHIFTIND = 0;
    static constexpr size_t size = 4;
    static constexpr size_t numberToSkipAsTuple = 1;

  private:
    device::array<Field<T, NDim>, 3> fs;
    const device::memory::host_string mName;
    LayoutStruct<NDim> mLayout;
  };
} // namespace TempLat

#endif
