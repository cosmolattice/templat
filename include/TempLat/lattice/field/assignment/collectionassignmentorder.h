#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNMENT_COLLECTIONASSIGNMENTORDER_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNMENT_COLLECTIONASSIGNMENTORDER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <concepts>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/getcomponent.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/lattice/latticeparameters.h"
#include "TempLat/lattice/memory/memorymanager.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"
#include "TempLat/util/log/saycomplete.h"
#include "TempLat/util/rangeiteration/for_in_range.h"

namespace TempLat
{
  namespace CollectionAssignmentDetail
  {
    template <typename X> using Naked = std::remove_cvref_t<X>;

    template <typename X>
    concept Composite = requires(const X &x) { x.assignedFields(); };
    template <typename X>
    concept PlainField = !Composite<X> && requires(const X &x) {
      x.getMemoryManager();
      x.getView();
    };
    template <typename X>
    concept Collection = !Composite<X> && !PlainField<X> && requires(const X &x) {
      Naked<X>::size;
      x.getComp(Tag<0>());
    };

    template <typename X> struct LeafOf {
      using type = void;
    };
    template <PlainField X> struct LeafOf<X> {
      using type = Naked<X>;
    };
    template <Composite X> struct LeafOf<X> {
      using type = Naked<decltype(std::declval<const X &>().assignedFields()[0])>;
    };

    /** @brief The field type an assignment to X writes, or void if unknown or not unique: X itself for a field, the
     *  component fields of a composite, those of the members of a collection. */
    template <typename X> constexpr auto leaf()
    {
      if constexpr (Collection<X>)
        return []<size_t... m>(std::index_sequence<m...>) {
          using First = decltype(leaf<decltype(std::declval<const X &>().getComp(Tag<0>()))>());
          if constexpr ((std::is_same_v<First,
                                        decltype(leaf<decltype(std::declval<const X &>().getComp(Tag<int(m)>()))>())> &&
                         ...))
            return First{};
          else
            return std::type_identity<void>{};
        }(std::make_index_sequence<Naked<X>::size>());
      else
        return std::type_identity<typename LeafOf<Naked<X>>::type>{};
    }
    template <typename X> using Leaf = typename decltype(leaf<X>())::type;

    /** @brief Appends the fields an assignment to x writes. */
    template <typename X, typename F> void appendWrittenFields(const X &x, std::vector<F> &out)
    {
      if constexpr (Composite<X>)
        for (const auto &f : x.assignedFields())
          out.push_back(f);
      else if constexpr (PlainField<X>)
        out.push_back(x);
      else
        for_in_range<0, Naked<X>::size>([&](auto m) { appendWrittenFields(x.getComp(m), out); });
    }

    /** @brief A new object of the type of x, for an assignment that cannot write x yet. */
    template <typename X> Naked<X> makeScratchLike(const X &x)
    {
      using F = Leaf<X>;
      std::vector<F> fields;
      appendWrittenFields(x, fields);
      using T = typename decltype(fields[0].getView())::value_type;
      const auto name = fields[0].getMemoryManager()->getName() + "_scratch";
      if constexpr (std::constructible_from<Naked<X>, std::string, decltype(fields[0].toolBox()),
                                            LatticeParameters<T>>) {
        static bool warned = false;
        if (!warned && MPICommReference::worldRank() == 0)
          sayShort << "TempLat warning: members of a collection read each other's targets in a cycle (e.g. a "
                      "swap). They are assigned through scratch fields, here "
                   << name << ".\n";
        warned = true;
        return Naked<X>(name, fields[0].toolBox(), LatticeParameters<T>());
      } else
        throw std::logic_error("Members of a collection read each other's targets in a cycle, which needs a scratch "
                               "copy of " +
                               name + ", but its type cannot be constructed from a name and a toolbox.");
    }

    template <typename R>
    concept List = requires {
      typename Naked<R>::Getter;
      Naked<R>::size;
    } && std::is_same_v<typename Naked<R>::Getter, GetComponent>;

    /** @brief Walks the expression r, or each member of the list r, bringing it to configuration space and updating
     *  its ghosts; reads of assignment targets get recorded. */
    template <typename R, typename Layout> void walk(R &r, const Layout &layout)
    {
      if constexpr (List<R>)
        for_in_range<0, Naked<R>::size>([&](auto k) {
          auto member = GetComponent::get(r, k);
          walk(member, layout);
        });
      else {
        ConfirmSpace::apply(r, layout, SpaceStateType::Configuration);
        GhostsHunter::apply(r);
      }
    }

    /** @brief The field type all members target(Tag<j>()), j < N, write, or void. */
    template <int N, typename Target> constexpr auto commonLeaf()
    {
      return []<size_t... m>(std::index_sequence<m...>) {
        using First = Leaf<decltype(std::declval<Target &>()(Tag<0>()))>;
        if constexpr ((std::is_same_v<First, Leaf<decltype(std::declval<Target &>()(Tag<int(m)>()))>> && ...))
          return std::type_identity<First>{};
        else
          return std::type_identity<void>{};
      }(std::make_index_sequence<N>());
    }

    /** @brief Calls f(Tag<j>()) for the runtime index j < N. */
    template <int N, typename Fun> void withTag(size_t j, Fun &&f)
    {
      for_in_range<0, N>([&](auto t) {
        if (decltype(t)::value == j) f(t);
      });
    }
  } // namespace CollectionAssignmentDetail

  /** @brief Assigns member j of a list to member j of a collection, j = 0 ... N - 1, one member at a time, such that no
   *  member's right-hand side reads a target already overwritten. target(Tag<j>()) gives member j, rhs(Tag<j>())
   *  the right-hand side, scratch(Tag<j>()) an object of the same type as member j.
   *
   *  Before assigning, each right-hand side is walked with all members' fields marked as targets. Member g reading a
   *  target of member h is assigned before h. Members whose reads form a cycle (e.g. a swap) cannot be ordered: one
   *  of them is assigned into scratch and copied back at the end. Without reads across members, they are assigned
   *  in order; the only cost is the walk. Member types whose written fields are unknown are assigned in
   *  order, unchecked.
   *
   *  Unit test: ctest -R test-collectionassignmentorder
   **/
  template <int N, typename Target, typename Rhs, typename Scratch>
  void assignMembersInOrder(Target &&target, Rhs &&rhs, Scratch &&scratch)
  {
    using namespace CollectionAssignmentDetail;
    // A named right-hand side: assigning a temporary collection would pick its move assignment, writing no data.
    auto assign = [&](auto &&to, auto j) {
      auto r = rhs(j);
      to = r;
    };
    auto assignInOrder = [&]() { for_in_range<0, N>([&](auto j) { assign(target(j), j); }); };

    if constexpr (N < 2)
      assignInOrder();
    else if constexpr (std::is_void_v<typename decltype(commonLeaf<N, Target>())::type>)
      assignInOrder();
    else {
      using F = typename decltype(commonLeaf<N, Target>())::type;
      // The fields of all members, and which member writes each.
      std::vector<F> fields;
      std::vector<size_t> owner;
      for_in_range<0, N>([&](auto j) {
        appendWrittenFields(target(j), fields);
        owner.resize(fields.size(), decltype(j)::value);
      });
      using Manager = Naked<decltype(*fields[0].getMemoryManager())>;
      std::vector<Manager *> managers;
      for (const auto &f : fields)
        managers.push_back(f.getMemoryManager().get());

      // reads[g * N + h]: the right-hand side of member g reads a target of member h != g.
      std::vector<char> reads(N * N, 0);
      bool readsOtherMember = false;
      for_in_range<0, N>([&](auto g) {
        auto r = rhs(g);
        typename Manager::AssignmentTargetGuard guard(managers);
        walk(r, fields[0].getLayout());
        for (size_t i = 0; i < fields.size(); ++i)
          if (owner[i] != decltype(g)::value && guard.wasRead(i)) {
            reads[decltype(g)::value * N + owner[i]] = 1;
            readsOtherMember = true;
          }
      });
      if (!readsOtherMember) {
        assignInOrder();
        return;
      }

      // Repeatedly assign the first member no remaining member reads; if there is none, defer the first remaining.
      std::vector<char> done(N, 0), deferred(N, 0);
      std::vector<size_t> order;
      for (size_t step = 0; step < size_t(N); ++step) {
        size_t next = N;
        for (size_t h = 0; h < size_t(N) && next == size_t(N); ++h) {
          bool isRead = false;
          for (size_t g = 0; g < size_t(N); ++g)
            isRead = isRead || (!done[g] && g != h && reads[g * N + h]);
          if (!done[h] && !isRead) next = h;
        }
        if (next == size_t(N)) {
          for (next = 0; done[next]; ++next) {
          }
          deferred[next] = 1;
        }
        done[next] = 1;
        order.push_back(next);
      }

      for (size_t j : order)
        withTag<N>(j, [&](auto t) {
          if (deferred[j])
            assign(scratch(t), t);
          else
            assign(target(t), t);
        });
      for (size_t j = 0; j < size_t(N); ++j)
        if (deferred[j])
          withTag<N>(j, [&](auto t) {
            std::vector<F> from, to;
            appendWrittenFields(scratch(t), from);
            appendWrittenFields(target(t), to);
            for (size_t k = 0; k < to.size(); ++k)
              to[k] = from[k];
          });
    }
  }
} // namespace TempLat

#endif
