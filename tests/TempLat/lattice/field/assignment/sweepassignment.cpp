/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/subtract.h"
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/util/constexpr_for.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/algebra/listoperators/listoperators.h"
#include "TempLat/lattice/algebra/listoperators/norm.h"

namespace TempLat
{
  /** @brief Assignments whose right-hand side reads the target at neighbouring sites give the same result as
   *  assigning into a separate field.
   **/
  struct SweepAssignmentTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace SweepAssignmentDetail
  {
    /** @brief The list whose member m is member m + 1 (mod N) of the collection c, shifted by S along axis 0. */
    template <int N, int S, typename C> struct Rotated {
      C c;
      static constexpr size_t size = N;
      using Getter = GetComponent;
      template <int M> auto getComp(Tag<M>) const
      {
        if constexpr (S == 0)
          return c.getComp(Tag<(M + 1) % N>());
        else
          return shift<S>(c.getComp(Tag<(M + 1) % N>()));
      }
    };
    template <int N, int S, typename C> Rotated<N, S, C> rotated(const C &c) { return {c}; }

    /** @brief Checks phi = expr(phi) in place against ref = expr(copy of phi), for the given band width
     *  (0: automatic). */
    template <size_t NDim, typename Expr>
    bool matchesOutOfPlace(const device::memory::host_ptr<MemoryToolBox<NDim>> &toolBox, device::Idx bandWidth,
                           Expr &&expr)
    {
      SpatialCoordinate<NDim> x(toolBox);
      Field<double, NDim> phi("phi", toolBox);
      Field<double, NDim> copy("copy", toolBox);
      Field<double, NDim> ref("ref", toolBox);

      // No symmetry, so that a read at the wrong site shows up.
      if constexpr (NDim == 3)
        phi = x(1_c) + 10. * x(2_c) + 100. * x(3_c) + x(1_c) * x(2_c) * x(3_c);
      else
        phi = x(1_c) + 10. * x(2_c) + x(1_c) * x(1_c) * x(2_c);
      copy = phi;

      toolBox->mSweepBuffer.setBandWidth(bandWidth);
      ref = expr(copy);
      phi = expr(phi);
      toolBox->mSweepBuffer.setBandWidth(0);

      return average((phi - ref) * (phi - ref)) == 0;
    }

    /** @brief As matchesOutOfPlace, for a composite field F with components first ... last. */
    template <typename F, int first, int last, typename Expr>
    bool compositeMatchesOutOfPlace(const device::memory::host_ptr<MemoryToolBox<3>> &toolBox, device::Idx bandWidth,
                                    Expr &&expr)
    {
      SpatialCoordinate<3> x(toolBox);
      F f("f", toolBox, LatticeParameters<double>());
      F copy("copy", toolBox, LatticeParameters<double>());
      F ref("ref", toolBox, LatticeParameters<double>());

      constexpr_for<first, last + 1>([&](auto _c) {
        constexpr Tag<(int)decltype(_c)::value> c;
        const double k = decltype(_c)::value + 1.;
        f(c) = k * x(1_c) + 10. * x(2_c) + 100. * x(3_c) + x(1_c) * x(2_c) * x(3_c) / k;
        copy(c) = f(c);
      });

      toolBox->mSweepBuffer.setBandWidth(bandWidth);
      ref = expr(copy);
      f = expr(f);
      toolBox->mSweepBuffer.setBandWidth(0);

      bool matches = true;
      constexpr_for<first, last + 1>([&](auto _c) {
        constexpr Tag<(int)decltype(_c)::value> c;
        matches = matches && average((f(c) - ref(c)) * (f(c) - ref(c))) == 0;
      });
      return matches;
    }

    /** @brief As matchesOutOfPlace, for a collection of NComp fields. */
    template <int NComp, typename Expr>
    bool collectionMatchesOutOfPlace(const device::memory::host_ptr<MemoryToolBox<3>> &toolBox, device::Idx bandWidth,
                                     Expr &&expr)
    {
      SpatialCoordinate<3> x(toolBox);
      FieldCollection<Field<double, 3>, NComp> f("f", toolBox), copy("copy", toolBox), ref("ref", toolBox);

      constexpr_for<0, NComp>([&](auto c) {
        const double k = decltype(c)::value + 1.;
        f(c) = k * x(1_c) + 10. * x(2_c) + 100. * x(3_c) + x(1_c) * x(2_c) * x(3_c) / k;
        copy(c) = f(c);
      });

      toolBox->mSweepBuffer.setBandWidth(bandWidth);
      ref = expr(copy);
      f = expr(f);
      toolBox->mSweepBuffer.setBandWidth(0);

      bool matches = true;
      constexpr_for<0, NComp>([&](auto c) { matches = matches && average((f(c) - ref(c)) * (f(c) - ref(c))) == 0; });
      return matches;
    }

    /** @brief As matchesOutOfPlace, for a collection of NMembers composite fields F with components first ... last. */
    template <typename F, int NMembers, int first, int last, typename Expr>
    bool compositeCollectionMatchesOutOfPlace(const device::memory::host_ptr<MemoryToolBox<3>> &toolBox,
                                              device::Idx bandWidth, Expr &&expr)
    {
      SpatialCoordinate<3> x(toolBox);
      FieldCollection<F, NMembers> f("f", toolBox), copy("copy", toolBox), ref("ref", toolBox);

      constexpr_for<0, NMembers>([&](auto m) {
        constexpr_for<first, last + 1>([&](auto _c) {
          constexpr Tag<(int)decltype(_c)::value> c;
          const double k = 10. * decltype(m)::value + decltype(_c)::value + 1.;
          f(m)(c) = k * x(1_c) + 10. * x(2_c) + 100. * x(3_c) + x(1_c) * x(2_c) * x(3_c) / k;
          copy(m)(c) = f(m)(c);
        });
      });

      toolBox->mSweepBuffer.setBandWidth(bandWidth);
      ref = expr(copy);
      f = expr(f);
      toolBox->mSweepBuffer.setBandWidth(0);

      bool matches = true;
      constexpr_for<0, NMembers>([&](auto m) {
        constexpr_for<first, last + 1>([&](auto _c) {
          constexpr Tag<(int)decltype(_c)::value> c;
          matches = matches && average((f(m)(c) - ref(m)(c)) * (f(m)(c) - ref(m)(c))) == 0;
        });
      });
      return matches;
    }
  } // namespace SweepAssignmentDetail

  inline void SweepAssignmentTester::Test(TDDAssertion &tdd)
  {
    using SweepAssignmentDetail::matchesOutOfPlace;
    constexpr size_t NDim = 3;
    auto toolBox = MemoryToolBox<NDim>::makeShared(8, 2);

    auto lapl = [](const auto &f) { return 2. * LatLapl(f); };
    auto lapl2 = [](const auto &f) { return f + 0.1 * LatLapl(LatLapl(f)); };
    auto down = [](const auto &f) { return shift<-2>(f); };
    auto up = [](const auto &f) { return shift<1>(f); };
    auto inPlane = [](const auto &f) { return shift<2>(f) + shift<-3>(f); };
    auto diagonal = [](const auto &f) { return shift<-1, 1, 0>(f) - f; };

    // Automatic band width.
    tdd.verify(matchesOutOfPlace(toolBox, 0, lapl));
    tdd.verify(matchesOutOfPlace(toolBox, 0, lapl2));
    tdd.verify(matchesOutOfPlace(toolBox, 0, down));
    tdd.verify(matchesOutOfPlace(toolBox, 0, up));
    tdd.verify(matchesOutOfPlace(toolBox, 0, inPlane));
    tdd.verify(matchesOutOfPlace(toolBox, 0, diagonal));

    // Band widths: one plane, a divisor of the local extent, a non-divisor (partial last band), and wider than
    // the local extent (a single step).
    for (device::Idx b : {1, 2, 3, 100}) {
      tdd.verify(matchesOutOfPlace(toolBox, b, lapl));
      tdd.verify(matchesOutOfPlace(toolBox, b, lapl2));
      tdd.verify(matchesOutOfPlace(toolBox, b, down));
      tdd.verify(matchesOutOfPlace(toolBox, b, up));
      tdd.verify(matchesOutOfPlace(toolBox, b, diagonal));
    }

    // Two dimensions.
    {
      auto toolBox2 = MemoryToolBox<2>::makeShared(8, 2);
      tdd.verify(matchesOutOfPlace(toolBox2, 0, lapl2));
      tdd.verify(matchesOutOfPlace(toolBox2, 3, down));
    }

    // Composite fields, assigned in a single kernel.
    {
      using SweepAssignmentDetail::compositeMatchesOutOfPlace;
      using C = ComplexField<double, NDim>;
      using U = SU2Field<double, NDim>;
      using D = SU2Doublet<double, NDim>;
      using A = SU2LieAlgebraField<double, NDim>;
      using H = SymTracelessField<double, NDim>;
      // A collection of composites is stored member by member, keeping one member's intermediates live at a time.
      static_assert(EvaluatesMembers<ComponentList<device::array<double, 4>, U, U>, NDim>);
      auto complexExpr = [](const auto &f) { return shift<-1>(f) * shift<2>(f) + f; };
      auto su2Expr = [](const auto &f) { return shift<-1>(f) * dagger(shift<2>(f)); };
      auto sum = [](const auto &f) { return shift<-1>(f) + shift<1>(f); };
      auto symExpr = [](const auto &f) { return shift<-1>(f) + 2. * shift<2>(f); };
      for (device::Idx b : {0, 3}) {
        tdd.verify(compositeMatchesOutOfPlace<C, 0, 1>(toolBox, b, complexExpr));
        tdd.verify(compositeMatchesOutOfPlace<U, 0, 3>(toolBox, b, su2Expr));
        tdd.verify(compositeMatchesOutOfPlace<D, 0, 3>(toolBox, b, sum));
        tdd.verify(compositeMatchesOutOfPlace<A, 1, 3>(toolBox, b, sum));
        tdd.verify(compositeMatchesOutOfPlace<H, 0, 4>(toolBox, b, symExpr));
      }
    }

    // Collections, assigned in a single kernel: components read each other in place and at neighbouring sites.
    {
      using SweepAssignmentDetail::collectionMatchesOutOfPlace;
      auto pointCross = [](auto &f) { return f * (1. + 1e-6 * norm2(f)); };
      auto shiftedCross = [](auto &f) { return f + 1e-6 * shift<-1>(norm2(f)) - LatLapl(f); };
      for (device::Idx b : {0, 3}) {
        tdd.verify(collectionMatchesOutOfPlace<2>(toolBox, b, shiftedCross));
        tdd.verify(collectionMatchesOutOfPlace<3>(toolBox, b, pointCross));
        tdd.verify(collectionMatchesOutOfPlace<3>(toolBox, b, shiftedCross));
      }
    }

    // Collections of composite fields, assigned in a single kernel.
    {
      using SweepAssignmentDetail::compositeCollectionMatchesOutOfPlace;
      auto sum = [](auto &f) { return shift<-1>(f) + shift<1>(f); };
      for (device::Idx b : {0, 3}) {
        tdd.verify(compositeCollectionMatchesOutOfPlace<SU2Field<double, NDim>, 2, 0, 3>(toolBox, b, sum));
        tdd.verify(compositeCollectionMatchesOutOfPlace<SU2LieAlgebraField<double, NDim>, 3, 1, 3>(toolBox, b, sum));
      }
    }

    // Members reading each other, in place and at neighbouring sites.
    {
      using SweepAssignmentDetail::collectionMatchesOutOfPlace;
      using SweepAssignmentDetail::compositeCollectionMatchesOutOfPlace;
      using SweepAssignmentDetail::rotated;
      auto inPlace = [](auto &f) { return rotated<3, 0>(f); };
      auto shifted = [](auto &f) { return rotated<3, -1>(f); };
      for (device::Idx b : {0, 3}) {
        tdd.verify(collectionMatchesOutOfPlace<3>(toolBox, b, inPlace));
        tdd.verify(collectionMatchesOutOfPlace<3>(toolBox, b, shifted));
        tdd.verify(compositeCollectionMatchesOutOfPlace<SU2Field<double, NDim>, 3, 0, 3>(toolBox, b, inPlace));
        tdd.verify(compositeCollectionMatchesOutOfPlace<SU2Field<double, NDim>, 3, 0, 3>(toolBox, b, shifted));
        tdd.verify(
            compositeCollectionMatchesOutOfPlace<SU2LieAlgebraField<double, NDim>, 3, 1, 3>(toolBox, b, inPlace));
        tdd.verify(
            compositeCollectionMatchesOutOfPlace<SU2LieAlgebraField<double, NDim>, 3, 1, 3>(toolBox, b, shifted));
      }
    }

    // The buffer is allocated only once an assignment needs it.
    {
      auto fresh = MemoryToolBox<NDim>::makeShared(8, 1);
      Field<double, NDim> phi("phi", fresh);
      Field<double, NDim> chi("chi", fresh);
      phi = 1.;
      chi = LatLapl(phi);
      phi = 2. * phi + chi;
      tdd.verify(!fresh->mSweepBuffer.isAllocated());
      phi = LatLapl(phi);
      tdd.verify(fresh->mSweepBuffer.isAllocated());
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SweepAssignmentTester> test;
}
