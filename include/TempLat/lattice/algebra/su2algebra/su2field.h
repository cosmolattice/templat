#ifndef COSMOINTERFACE_SU2ALGEBRA_SU2FIELD_H
#define COSMOINTERFACE_SU2ALGEBRA_SU2FIELD_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/field/assignablefieldcollection.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2get.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"

#include "TempLat/parallel/device.h"

#include "TempLat/lattice/algebra/su2algebra/su2binaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2commutator.h"

namespace TempLat
{
  /** @brief A class which implements a SU2 field (group).
   * We choose a representation with 4 real numbers per lattice site, c0, c1, c2, c3.
   * These are subject to the unitarity constraint c0^2 + c1^2 + c2^2 + c3^2 = 1, which can be enforced by the
   * unitarize() method.
   * A single link variable is then represented as U = c0 * I + i * (c1 * sigma1 + c2 * sigma2 + c3 * sigma3), where
   * sigma1, sigma2, sigma3 are the Pauli matrices.
   *
   * Unit test: ctest -R test-su2field
   *
   * @vocab-summary A single $SU(2)$ link $U_i$ for fixed $i$, held as four reals $c_0\dots c_3$ with $U =
   * c_0\mathbb{1} + i\,c_a\sigma_a$ and $\sum_a c_a^2 = 1$. unitarize() restores the constraint after
   * accumulated rounding.
   * @vocab-signature SU2Field<T, NDim> U("U", toolBox);
   **/
  template <typename T, size_t _NDim = 0> class SU2Field
  {
  public:
    // Put public methods here. These should change very little over time.
    static_assert(_NDim != 0, "NDim template parameter is required. Use e.g. SU2Field<double, 3>.");

    static constexpr size_t NDim = _NDim;

    SU2Field(Field<T, NDim> f0, Field<T, NDim> f1, Field<T, NDim> f2, Field<T, NDim> f3)
        : fs{{f0, f1, f2, f3}}, mName("NoName"), mLayout(fs[0].getToolBox()->mLayouts.getConfigSpaceLayout())
    {
    }

    SU2Field(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
             LatticeParameters<T> pLatPar = LatticeParameters<T>())
        : fs{{
              Field<T, NDim>(name + "_0", toolBox, pLatPar), //
              Field<T, NDim>(name + "_1", toolBox, pLatPar), //
              Field<T, NDim>(name + "_2", toolBox, pLatPar), //
              Field<T, NDim>(name + "_3", toolBox, pLatPar)  //
          }},
          mName(name), mLayout(toolBox->mLayouts.getConfigSpaceLayout())
    {
      fs[0] = T(1);
      fs[0].updateGhosts();
    }

    template <int N> const auto &SU2Get(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 3, "SU2Get: N must be between 0 and 3 for SU2Field");
      return fs[N];
    }

    template <int M> auto &operator()(Tag<M> t)
    {
      static_assert(M >= 0 && M <= 3, "Operator(): M must be between 0 and 3 for SU2Field");
      return fs[M];
    }

    template <int M> const auto &operator()(Tag<M> t) const
    {
      static_assert(M >= 0 && M <= 3, "Operator(): M must be between 0 and 3 for SU2Field");
      return fs[M];
    }

    template <typename R> void operator=(R &&r)
    {
      // Confirm config-space / ghost requirements by walking the WHOLE SU(2) expression r (linear in its
      // structure), instead of building r.SU2Get(k) per component (the exponential per-component
      // expansion, instantiated only to be walked then discarded). One call per target component is
      // kept: onBeforeAssignment also confirms the target field fs[k] and flags its host mirror.
      fs[0].onBeforeAssignment(r);
      fs[1].onBeforeAssignment(r);
      fs[2].onBeforeAssignment(r);
      fs[3].onBeforeAssignment(r);

      PreGet::apply(r);

      const auto view0 = fs[0].getView();
      const auto view1 = fs[1].getView();
      const auto view2 = fs[2].getView();
      const auto view3 = fs[3].getView();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply(
            [&](const auto &...args) {
              auto result = DoEval::eval(r, args...);
              view0(args...) = result[0];
              view1(args...) = result[1];
              view2(args...) = result[2];
              view3(args...) = result[3];
            },
            idx);
      };
      device::iteration::foreach ("SU2ConfigViewAssign", mLayout, functor);

      PostGet::apply(r);

      fs[0].setGhostsAreStale();
      fs[1].setGhostsAreStale();
      fs[2].setGhostsAreStale();
      fs[3].setGhostsAreStale();
    }

    // Recompute c0 from c1,c2,c3 to enforce SU(2) unitarity constraint.
    void unitarize()
    {
      const auto view1 = fs[1].getView();
      const auto view2 = fs[2].getView();
      const auto view3 = fs[3].getView();
      const auto view0 = fs[0].getView();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::apply(
            [&](const auto &...args) {
              T c1 = view1(args...);
              T c2 = view2(args...);
              T c3 = view3(args...);
              view0(args...) = sqrt(T(1) - c1 * c1 - c2 * c2 - c3 * c3);
            },
            idx);
      };
      device::iteration::foreach ("SU2Unitarize", mLayout, functor);

      fs[0].setGhostsAreStale();
      fs[0].updateGhosts();
    }

    std::string toString() const { return *mName; }

    auto getDx() const { return GetDx::getDx(fs[0]); }
    auto getKIR() const { return GetKIR::getKIR(fs[0]); }

    inline auto getToolBox() const { return GetToolBox::get(fs[0]); }

    // Space/ghost confirmation forwarded to the 4 component fields, so this leaf can be walked as a whole
    // SU(2) expression (see operator=). GhostsHunter::apply on a Field is a no-op (matches the prior
    // per-component path where r.SU2Get(k) was the bare field fs[k]).
    void doWeNeedGhosts() const
    {
      GhostsHunter::apply(fs[0]);
      GhostsHunter::apply(fs[1]);
      GhostsHunter::apply(fs[2]);
      GhostsHunter::apply(fs[3]);
    }
    device::Idx confirmGhostsUpToDate() const
    {
      MemoryManager<T, NDim> *mgrs[] = {fs[0].getMemoryManager().get(), fs[1].getMemoryManager().get(),
                                        fs[2].getMemoryManager().get(), fs[3].getMemoryManager().get()};
      return MemoryManager<T, NDim>::confirmGhostsUpToDateBatch(mgrs);
    }
    void confirmSpace(const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      ConfirmSpace::apply(fs[0], newLayout, spaceType);
      ConfirmSpace::apply(fs[1], newLayout, spaceType);
      ConfirmSpace::apply(fs[2], newLayout, spaceType);
      ConfirmSpace::apply(fs[3], newLayout, spaceType);
    }

    inline void updateGhosts()
    {
      MemoryManager<T, NDim> *mgrs[] = {fs[0].getMemoryManager().get(), fs[1].getMemoryManager().get(),
                                        fs[2].getMemoryManager().get(), fs[3].getMemoryManager().get()};
      MemoryManager<T, NDim>::updateGhostsBatch(mgrs);
    }

    template <typename... IDX>
      requires requires(Field<T, NDim> f, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(f, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<T, 4> result;
      result[0] = fs[0].eval(idx...);
      result[1] = fs[1].eval(idx...);
      result[2] = fs[2].eval(idx...);
      result[3] = fs[3].eval(idx...);
      return result;
    }

    using Getter = SU2Getter;

    static constexpr size_t SHIFTIND = 0;
    static constexpr size_t size = 4;
    static constexpr size_t numberToSkipAsTuple = 0;

  protected:
    device::array<Field<T, NDim>, 4> fs;
    const device::memory::host_string mName;
    LayoutStruct<NDim> mLayout;
  };
} // namespace TempLat

#endif
