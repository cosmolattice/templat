#ifndef COSMOINTERFACE_SU2ALGEBRA_SU2DOUBLET_H
#define COSMOINTERFACE_SU2ALGEBRA_SU2DOUBLET_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getdx.h"
#include "TempLat/lattice/algebra/helpers/getkir.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/field/assignablefieldcollection.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2doubletget.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"

#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief A class which implements su2doublets.
   *
   *
   * Unit test: ctest -R test-su2doublet
   **/
  template <typename T, size_t _NDim = 0> class SU2Doublet
  {
  public:
    // Put public methods here. These should change very little over time.
    static_assert(_NDim != 0, "NDim template parameter is required. Use e.g. SU2Doublet<double, 3>.");

    static constexpr size_t NDim = _NDim;

    SU2Doublet(Field<T, NDim> f1, Field<T, NDim> f2, Field<T, NDim> f3, Field<T, NDim> f4)
        : fs{{f1, f2, f3, f4}}, mName("NoName"), mLayout(f1.getToolBox()->mLayouts.getConfigSpaceLayout())
    {
    }
    SU2Doublet(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
               LatticeParameters<T> pLatPar = LatticeParameters<T>())
        : mName(name), fs{{
                           Field<T, NDim>(name + "_0", toolBox, pLatPar), //
                           Field<T, NDim>(name + "_1", toolBox, pLatPar), //
                           Field<T, NDim>(name + "_2", toolBox, pLatPar), //
                           Field<T, NDim>(name + "_3", toolBox, pLatPar)  //
                       }},
          mLayout(toolBox->mLayouts.getConfigSpaceLayout())
    {
    }

    template <int N> const Field<T, NDim> &SU2DoubletGet(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "SU2DoubletGet: N must be between 0 and 3 for SU2Doublet");
      return fs[t];
    }

    template <typename... IDX>
      requires requires(Field<T, NDim> f, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(f, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<T, 4> result;
      result[0] = DoEval::eval(fs[0], idx...);
      result[1] = DoEval::eval(fs[1], idx...);
      result[2] = DoEval::eval(fs[2], idx...);
      result[3] = DoEval::eval(fs[3], idx...);
      return result;
    }

    template <int M> auto &operator()(Tag<M> t)
    {
      static_assert(M >= 0 && M <= 3, "Operator(): M must be between 0 and 3 for SU2Doublet");
      return fs[t];
    }

    template <typename R> void operator=(R &&r)
    {
      fs[0].onBeforeAssignment(std::decay_t<R>::Getter::get(r, 0_c));
      fs[1].onBeforeAssignment(std::decay_t<R>::Getter::get(r, 1_c));
      fs[2].onBeforeAssignment(std::decay_t<R>::Getter::get(r, 2_c));
      fs[3].onBeforeAssignment(std::decay_t<R>::Getter::get(r, 3_c));

      PreGet::apply(r);

      const auto views = device::make_tuple(fs[0].getView(), fs[1].getView(), fs[2].getView(), fs[3].getView());

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply(
            [&](auto &&...args) {
              auto result = DoEval::eval(r, args...);
              constexpr_for<0, size>([&](auto _i) {
                constexpr size_t i = decltype(_i)::value;
                device::get<i>(views)(args...) = result[i];
              });
            },
            idx);
      };
      device::iteration::foreach ("SU2DoubleConfigViewAssign", mLayout, functor);

      PostGet::apply(r);

      fs[0].setGhostsAreStale();
      fs[1].setGhostsAreStale();
      fs[2].setGhostsAreStale();
      fs[3].setGhostsAreStale();
    }

    template <typename R> void operator+=(R &&r) { (*this) = (*this) + r; }

    std::string toString() const { return *mName; }

    device::memory::host_ptr<MemoryToolBox<NDim>> getToolBox() const { return GetToolBox::get(fs[0]); }

    auto getDx() const { return GetDx::getDx(fs[0]); }
    auto getKIR() const { return GetKIR::getKIR(fs[0]); }

    inline void updateGhosts()
    {
      fs[0].updateGhosts();
      fs[1].updateGhosts();
      fs[2].updateGhosts();
      fs[3].updateGhosts();
    }

    using Getter = SU2DoubletGetter;
    static constexpr size_t size = 4;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    const device::memory::host_string mName;

    device::array<Field<T, NDim>, 4> fs;

    LayoutStruct<NDim> mLayout;
  };
} // namespace TempLat

#endif
