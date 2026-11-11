#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNABLEFIELDCOLLECTION_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNABLEFIELDCOLLECTION_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include <optional>

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/field/assignment/collectionassignmentorder.h"
#include "TempLat/util/shiftedindexsequence.h"
#include "TempLat/util/foreach.h"
#include "TempLat/util/istuplelike.h"
#include "TempLat/util/tuple_size.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/tuplemaker.h"
#include "TempLat/util/flattentuple.h"

namespace TempLat
{
  // if want to try to have dynamical access somehow, see
  // https://www.justsoftwaresolutions.co.uk/cplusplus/getting-tuple-elements-with-runtime-index.html
  // seems hardly compatible with indexed fields though.

  /** @brief A class which implements collections of object which are assignable.
   *
   *
   * Unit test: ctest -R test-assignablefieldcollection
   **/
  template <class Q, class... Args> class AssignableCollectionBase
  {
  public:
    static constexpr size_t size = sizeof...(Args);

    AssignableCollectionBase(Args... args)
        : // This constructor is used to create composite structure such as SU(N) matrices, from fields.
          fs(std::make_tuple(args...))
    {
    }

    void operator=(const AssignableCollectionBase<Q, Args...> &other)
    {
      for_in_range<0, size>([&](auto j) { std::get<j>(fs) = std::get<j>(other.fs); });
    }

    /** @brief Assigns the list r member by member, in an order in which no member reads another member's target
     *  after it was overwritten, see assignMembersInOrder. */
    template <typename R>
      requires(!std::is_same_v<R, AssignableCollectionBase<Q, Args...>>)
    void operator=(R &&r)
    {
      using nakedR = std::decay_t<R>;
      // Scratch members are allocated only for a cycle of reads; tuples assigned to are usually temporaries.
      std::tuple<std::optional<Args>...> scratch;
      assignMembersInOrder<size>([&](auto j) -> auto & { return std::get<j>(fs); },
                                 [&](auto j) { return nakedR::Getter::get(r, j); },
                                 [&](auto j) -> auto & {
                                   auto &s = std::get<j>(scratch);
                                   if (!s) s.emplace(CollectionAssignmentDetail::makeScratchLike(std::get<j>(fs)));
                                   return *s;
                                 });
    }

    template <int M> auto &operator()(Tag<M> t)
    {
      static_assert(M >= Q::SHIFTIND && M < size + Q::SHIFTIND, "Index out of bounds in field collection.");
      // Also work with operator [], but might be confusing that it is NOT an int (cant call f[1] -> f[1_c])
      return std::get<M - Q::SHIFTIND>(fs);
    }

    template <int M> const auto &operator()(Tag<M> t) const
    {
      static_assert(M >= Q::SHIFTIND && M < size + Q::SHIFTIND, "Index out of bounds in field collection.");
      return std::get<M - Q::SHIFTIND>(fs);
    }

  protected:
    std::tuple<Args...> fs;
  };

  template <class Q> // Specialise to empty field collection, that does nothing. Useful to define general models for
                     // example, which may have some empty collections.
  class AssignableCollectionBase<Q>
  {
  public:
    // Put public methods here. These should change very little over time.
    AssignableCollectionBase() = default;
    template <int N> int operator()(Tag<N> t) { return 0; }
    template <int N> int getComp(Tag<N> t) { return 0; }
    template <typename R> void operator=(R &&g) {}
    static constexpr size_t size = 0;
  };
} // namespace TempLat

#endif
