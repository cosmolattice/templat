#ifndef TEMPLAT_LATTICE_FIELD_COLLECTIONS_FIELDCOLLECTION_H
#define TEMPLAT_LATTICE_FIELD_COLLECTIONS_FIELDCOLLECTION_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

#include "TempLat/util/latinindiceslist.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/field/assignment/assignviews.h"
#include "TempLat/lattice/field/assignment/componentlist.h"
#include "TempLat/lattice/field/assignment/prepareassignment.h"

#include <string>
#include <tuple>

namespace TempLat
{
  /** @brief A class which
   * Field collections. Allows to have collection of composite vector fields as well.
   *
   *
   * Unit test: ctest -R test-fieldcollection
   *
   * @vocab-summary $N$ fields of any one type under a single tag-indexed name.
   * @vocab-signature FieldCollection<FieldType, N> coll("coll", toolBox);
   * @vocab-tags Collection
   **/
  template <class Arg, int N, bool flatAssign = false, int SHIFTIND = 0> class FieldCollection
  {
  public:
    static constexpr size_t NDim = GetNDim::get<Arg>();
    using RT = GetGetReturnType<Arg>::type;
    using T = typename GetFloatType<RT>::type;

    FieldCollection(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
                    LatticeParameters<T> pLatPar = LatticeParameters<T>())
    {
      for (int i = 0; i < N; ++i) {
        fs.push_back(Arg(name + "_" + std::to_string(i + SHIFTIND), toolBox, pLatPar));
      }
    }

    FieldCollection(const FieldCollection &) = default;
    FieldCollection(FieldCollection &&) noexcept = default;
    FieldCollection &operator=(FieldCollection &&) noexcept = default;
    ~FieldCollection() = default;

    template <int M>
      requires(M >= SHIFTIND && M < N + SHIFTIND)
    auto operator()(Tag<M> t) const
    {
      return fs[t - Tag<SHIFTIND>()];
    }
    template <int M>
      requires(M >= SHIFTIND && M < N + SHIFTIND)
    auto operator[](Tag<M> t) const
    {
      return fs[t - Tag<SHIFTIND>()];
    }

    template <int M> auto getComp(Tag<M> t) const { return fs[t]; }

    /** @brief Assigns the list r. Collections of plain fields, and of composite fields exposing assignedFields, assign
     *  all members in one kernel, so that members may read each other, see assignViews; others assign member by
     *  member. */
    template <typename R> void operator=(R &&r)
    {
      using nakedR = std::remove_reference_t<R>;
      if constexpr (nakedR::size == N &&
                    (requires(const Arg &a) { a.getView(); } || requires(const Arg &a) { a.assignedFields(); }))
        assignAtOnce(r, std::make_index_sequence<N>());
      else
        for_in_range<0, nakedR::size>([&](auto i) { (*this).getComp(i) = nakedR::Getter::get(r, i); });
    }

    template <typename R> void operator+=(R &&r) { (*this) = (*this) + r; }

    void operator=(const FieldCollection<Arg, N, flatAssign, SHIFTIND> &other)
    { // overwrite the default = operator.
      for (size_t i = 0; i < fs.size(); ++i)
        fs[i] = other.fs[i];
    }

    std::string toString(device::Idx i) const { return fs[i - SHIFTIND].toString(); }

    using Getter = GetComponent;
    static constexpr size_t size = N;

  private:
    template <typename A> struct FieldCount;
    template <typename F, size_t K> struct FieldCount<device::array<F, K>> {
      static constexpr size_t value = K;
    };

    /** @brief The fields an assignment of f writes, and the value component the first one receives. */
    static auto assignedFieldsOf(const Arg &f)
    {
      if constexpr (requires { f.assignedFields(); })
        return f.assignedFields();
      else
        return device::array<Arg, 1>{f};
    }
    static constexpr size_t firstAssignedComponent()
    {
      if constexpr (requires { Arg::firstAssignedComponent; })
        return Arg::firstAssignedComponent;
      else
        return 0;
    }

    template <typename R, size_t... i> void assignAtOnce(R &r, std::index_sequence<i...>)
    {
      using nakedR = std::remove_reference_t<R>;
      ComponentList<T, decltype(nakedR::Getter::get(r, Tag<int(i)>()))...> g(nakedR::Getter::get(r, Tag<int(i)>())...);

      // Member m writes the K fields m * K + k.
      constexpr size_t K = FieldCount<std::remove_cvref_t<decltype(assignedFieldsOf(fs[0]))>>::value;
      const auto first = assignedFieldsOf(fs[0])[0];
      device::array<decltype(first.getView()), N * K> views;
      std::array<MemoryManager<T, NDim> *, N * K> managers;
      for (size_t m = 0; m < N; ++m)
        for (size_t k = 0; k < K; ++k) {
          const auto f = assignedFieldsOf(fs[m])[k];
          views[m * K + k] = f.getView();
          managers[m * K + k] = f.getMemoryManager().get();
        }

      const auto &layout = first.getLayout();
      GhostReach targetReach;
      bool readsOtherMember = false;
      if constexpr (K == 1)
        targetReach = prepareAssignment<T, NDim>(g, managers, layout);
      else
        std::tie(targetReach, readsOtherMember) = prepareCollectionAssignment<T, NDim>(g, managers, K, layout);
      assignViews<firstAssignedComponent()>("CollectionAssign", views, g, layout, targetReach,
                                            first.toolBox()->mSweepBuffer, readsOtherMember);
      for (auto *m : managers)
        m->setGhostsAreStale();
    }

    std::vector<Arg> fs;
  };
} // namespace TempLat

#endif
