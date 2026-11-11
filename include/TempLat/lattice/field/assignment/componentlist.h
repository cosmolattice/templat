#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNMENT_COMPONENTLIST_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNMENT_COMPONENTLIST_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/helpers/confirmghosts.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/parallel/device.h"

namespace TempLat
{
  /** @brief The member expressions of a list right-hand side as one expression, for a collection assigning all
   *  members in a single kernel (FieldCollection::operator=). Its value at a site is the array of the member
   *  values; evalMember gives one member's value, for kernels storing member by member.
   *
   *  Unit test: ctest -R test-sweepassignment
   **/
  template <typename T, typename... Rs> class ComponentList
  {
  public:
    explicit ComponentList(const Rs &...rs) : mRs(rs...) {}

    /** @brief Marks the assignment kernel as one evaluating many fields per site, see FusesComponents. */
    static constexpr bool fusesComponents = true;
    static constexpr size_t size = sizeof...(Rs);

    template <typename... IDX> DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      return device::apply(
          [&](const auto &...r) { return device::array<T, sizeof...(Rs)>{T(DoEval::eval(r, idx...))...}; }, mRs);
    }

    template <size_t m, typename... IDX> DEVICE_INLINE_FUNCTION auto evalMember(const IDX &...idx) const
    {
      return DoEval::eval(device::get<m>(mRs), idx...);
    }

    /** @brief Calls f on every member expression, in order. Host only. */
    template <typename F> void forEachMember(F &&f)
    {
      device::apply([&](auto &...r) { (f(r), ...); }, mRs);
    }

    device::Idx confirmGhostsUpToDate(const GhostReach &reach) const
    {
      return device::apply([&](const auto &...r) { return (device::Idx(0) + ... + ConfirmGhosts::apply(r, reach)); },
                           mRs);
    }

    template <size_t LayoutNDim>
    void confirmSpace(const LayoutStruct<LayoutNDim> &newLayout, const SpaceStateType &spaceType) const
    {
      device::apply([&](const auto &...r) { (ConfirmSpace::apply(r, newLayout, spaceType), ...); }, mRs);
    }

    void preGet()
    {
      device::apply([](auto &...r) { (PreGet::apply(r), ...); }, mRs);
    }

    void postGet()
    {
      device::apply([](auto &...r) { (PostGet::apply(r), ...); }, mRs);
    }

  private:
    device::tuple<Rs...> mRs;
  };
} // namespace TempLat

#endif
