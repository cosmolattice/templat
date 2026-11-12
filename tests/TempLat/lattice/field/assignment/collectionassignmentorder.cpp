/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/field/assignment/collectionassignmentorder.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/field/collections/vectorfieldcollection.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/listoperators/listoperators.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/util/assignabletuple.h"
#include "TempLat/util/constexpr_for.h"
#include "TempLat/util/rangeiteration/for_in_range.h"

namespace TempLat
{
  /** @brief Collections assigned member by member (VectorFieldCollection, AssignableTuple) give the same result as
   *  assigning into a separate collection, also when a member reads another member's target.
   **/
  struct CollectionAssignmentOrderTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace CollectionAssignmentOrderDetail
  {
    constexpr size_t NDim = 3;

    /** @brief The largest mean square difference between the component fields of the groups of a and b. */
    template <int G, int first, int last, typename A, typename B> double deviation(const A &a, const B &b)
    {
      double dev = 0;
      for_in_range<0, G>([&](auto g) {
        for_in_range<1, NDim + 1>([&](auto i) {
          if constexpr (first == last)
            dev = std::max(dev, (double)average(pow<2>(a(g)(i) - b(g)(i))));
          else
            constexpr_for<first, last + 1>([&](auto _c) {
              constexpr Tag<(int)decltype(_c)::value> c;
              dev = std::max(dev, (double)average(pow<2>(a(g)(i)(c)-b(g)(i)(c))));
            });
        });
      });
      return dev;
    }

    /** @brief Fills every component field of the groups of w with a different function without symmetry. */
    template <int G, int first, int last, typename W, typename X> void fill(W &w, const X &x)
    {
      for_in_range<0, G>([&](auto g) {
        for_in_range<1, NDim + 1>([&](auto i) {
          constexpr_for<first, last + 1>([&](auto _c) {
            const double k = 1. + decltype(g)::value + 0.3 * decltype(i)::value + 0.1 * decltype(_c)::value;
            auto value = cos(0.4 * k * x(1_c) - 0.3 * x(2_c) + 0.2 * k * x(3_c)) + 0.1 * k;
            if constexpr (first == last)
              w(g)(i) = value;
            else
              w(g)(i)(Tag<(int)decltype(_c)::value>()) = value;
          });
        });
      });
    }

    /** @brief Checks w = expr(w) in place against ref = expr(w) for a VectorFieldCollection of G groups. */
    template <typename F, int G, int first, int last>
    bool groupsMatchOutOfPlace(const device::memory::host_ptr<MemoryToolBox<NDim>> &toolBox, auto &&expr)
    {
      SpatialCoordinate<NDim> x(toolBox);
      VectorFieldCollection<F, G> w("w", toolBox), ref("ref", toolBox);
      fill<G, first, last>(w, x);
      ref = expr(w);
      w = expr(w);
      return deviation<G, first, last>(w, ref) < 1e-26;
    }
  } // namespace CollectionAssignmentOrderDetail

  inline void CollectionAssignmentOrderTester::Test(TDDAssertion &tdd)
  {
    using namespace CollectionAssignmentOrderDetail;
    auto toolBox = MemoryToolBox<NDim>::makeShared(8, 1);
    using Plain = Field<double, NDim>;
    using Lie = SU2LieAlgebraField<double, NDim>;

    // No member reads another: assigned in order.
    tdd.verify(groupsMatchOutOfPlace<Plain, 2, 0, 0>(
        toolBox, [](const auto &w) { return make_list(2. * w(0_c), w(1_c) + shift<1>(w(1_c))); }));
    // Group 1 reads group 0: assigned first.
    tdd.verify(groupsMatchOutOfPlace<Plain, 2, 0, 0>(
        toolBox, [](const auto &w) { return make_list(2. * w(0_c), w(0_c) + shift<1>(w(1_c))); }));
    // Swap, and a rotation of three groups with shifts: cycles, one group goes through scratch fields.
    tdd.verify(groupsMatchOutOfPlace<Plain, 2, 0, 0>(toolBox, [](const auto &w) { return make_list(w(1_c), w(0_c)); }));
    tdd.verify(groupsMatchOutOfPlace<Plain, 3, 0, 0>(
        toolBox, [](const auto &w) { return make_list(w(1_c), shift<-1>(w(2_c)), w(0_c) + 0.5 * w(1_c)); }));
    // Composite members.
    tdd.verify(groupsMatchOutOfPlace<Lie, 2, 1, 3>(
        toolBox, [](const auto &w) { return make_list(shift<-1>(w(1_c)), w(0_c) + w(1_c)); }));

    // Tuples of fields.
    {
      SpatialCoordinate<NDim> x(toolBox);
      Plain f("f", toolBox), g("g", toolBox), refF("refF", toolBox), refG("refG", toolBox);
      auto setUp = [&]() {
        f = cos(0.3 * x(1_c) + 0.2 * x(2_c));
        g = sin(0.5 * x(2_c) - 0.4 * x(3_c));
      };
      auto matches = [&]() { return average(pow<2>(f - refF)) + average(pow<2>(g - refG)) < 1e-26; };

      setUp();
      make_list(refF, refG) = make_list(g, shift<1>(f));
      make_list(f, g) = make_list(g, shift<1>(f));
      tdd.verify(matches());

      setUp();
      make_list(refF, refG) = make_list(2. * f, f + g);
      make_list(f, g) = make_list(2. * f, f + g);
      tdd.verify(matches());
    }
  }
} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::CollectionAssignmentOrderTester> test;
}
