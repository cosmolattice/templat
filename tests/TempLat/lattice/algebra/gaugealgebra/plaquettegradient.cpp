/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include "TempLat/lattice/algebra/gaugealgebra/plaquettegradient.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/gaugealgebra/u1field.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/listoperators/listoperators.h"
#include "TempLat/lattice/field/collections/vectorfield.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/util/constexpr_for.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"

namespace TempLat
{
  /** @brief plaqGrad equals sum_{j != i} [plaq - plaqBack], member by member and when a collection assignment
   *  evaluates it for all directions at once.
   **/
  struct PlaquetteGradientTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace PlaquetteGradientDetail
  {
    constexpr size_t NDim = 3;

    /** @brief The largest mean square difference between the components first ... last of the collections a, b. */
    template <int first, int last, typename A, typename B> double deviation(A &a, B &b)
    {
      double dev = 0;
      for_in_range<1, NDim + 1>([&](auto i) {
        constexpr_for<first, last + 1>([&](auto _c) {
          constexpr Tag<(int)decltype(_c)::value> c;
          dev = std::max(dev, (double)average((a(i)(c) - b(i)(c)) * (a(i)(c) - b(i)(c))));
        });
      });
      return dev;
    }

    template <typename U> auto explicitGrad(const U &Us, auto i)
    {
      return Total(j, 1, (int)NDim, IfElse(i != j, plaq(Us, i, j) - plaqBack(Us, i, j), ZeroType()));
    }
  } // namespace PlaquetteGradientDetail

  inline void PlaquetteGradientTester::Test(TDDAssertion &tdd)
  {
    using namespace PlaquetteGradientDetail;
    using T = double;
    auto toolBox = MemoryToolBox<NDim>::makeShared(8, 1);
    SpatialCoordinate<NDim> x(toolBox);

    // SU(2) links without symmetry; plaqGrad's identities hold for any quaternions.
    VectorField<SU2Field<T, NDim>> U("U", toolBox);
    for_in_range<1, NDim + 1>([&](auto i) {
      constexpr_for<0, 4>([&](auto _c) {
        constexpr Tag<(int)decltype(_c)::value> c;
        const double k = 4. * decltype(i)::value + decltype(_c)::value;
        U(i)(c) = cos(0.3 * k * x(1_c) + 0.7 * x(2_c) - 0.2 * k * x(3_c)) + 0.1 * k;
      });
    });

    // Member by member: the explicit expression.
    {
      VectorField<SU2Field<T, NDim>> member("member", toolBox), ref("ref", toolBox);
      for_in_range<1, NDim + 1>([&](auto i) {
        member(i) = GetComponent::get(plaqGrad(U), Tag<decltype(i)::value - 1>());
        ref(i) = explicitGrad(U, i);
      });
      tdd.verify(deviation<0, 3>(member, ref) == 0);
    }

    // Whole: a collection assignment computes each forward plaquette once.
    {
      VectorField<SU2Field<T, NDim>> whole("whole", toolBox), ref("ref", toolBox);
      static_assert(wholeListShares<decltype(plaqGrad(U))>());
      whole = plaqGrad(U);
      for_in_range<1, NDim + 1>([&](auto i) { ref(i) = explicitGrad(U, i); });
      tdd.verify(deviation<0, 3>(whole, ref) < 1e-26);
    }

    // The kick, and a mix with a list that is evaluated member by member.
    {
      VectorField<SU2LieAlgebraField<T, NDim>> Pi("Pi", toolBox), ref("ref", toolBox);
      for_in_range<1, NDim + 1>([&](auto i) {
        constexpr_for<1, 4>([&](auto _c) {
          constexpr Tag<(int)decltype(_c)::value> c;
          Pi(i)(c) = 0.1 * decltype(i)::value * x(1_c) - 0.05 * decltype(_c)::value * x(3_c);
          ref(i)(c) = Pi(i)(c);
        });
      });
      auto kick = Pi - 0.01 * plaqGrad(U);
      static_assert(decltype(wholeOf<T, NDim, NDim>(kick))::shares);
      static_assert(!EvaluatesMembers<WholeAssignment<decltype(wholeOf<T, NDim, NDim>(kick))>, NDim>);
      Pi = kick;
      for_in_range<1, NDim + 1>([&](auto i) { ref(i) = ref(i) - 0.01 * explicitGrad(U, i); });
      tdd.verify(deviation<1, 3>(Pi, ref) < 1e-26);

      // An identically zero summand keeps the sharing.
      auto withZero = Pi - 0.01 * plaqGrad(U) - 0.5 * MakeVector(i, 1, NDim, MakeSU2(a, ZeroType()));
      static_assert(decltype(wholeOf<T, NDim, NDim>(withZero))::shares);
      Pi = withZero;
      for_in_range<1, NDim + 1>([&](auto i) { ref(i) = ref(i) - 0.01 * explicitGrad(U, i); });
      tdd.verify(deviation<1, 3>(Pi, ref) < 1e-26);

      Pi = Pi + 0.3 * plaqGrad(U) - MakeVector(i, 1, NDim, 0.5 * Pi(i));
      for_in_range<1, NDim + 1>([&](auto i) { ref(i) = ref(i) + 0.3 * explicitGrad(U, i) - 0.5 * ref(i); });
      tdd.verify(deviation<1, 3>(Pi, ref) < 1e-26);
    }

    // U(1) links.
    {
      VectorField<U1Field<T, NDim>> V("V", toolBox);
      for_in_range<1, NDim + 1>(
          [&](auto i) { V(i) = complexPhase(0.4 * decltype(i)::value * x(1_c) + 0.3 * x(2_c) - 0.2 * x(3_c)); });
      VectorField<ComplexField<T, NDim>> whole("whole", toolBox), ref("ref", toolBox);
      whole = plaqGrad(V);
      for_in_range<1, NDim + 1>([&](auto i) { ref(i) = explicitGrad(V, i); });
      tdd.verify(deviation<0, 1>(whole, ref) < 1e-26);
    }
  }
} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::PlaquetteGradientTester> test;
}
