/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/algebra/helpers/ghostreach.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/lattice/algebra/spatialderivatives/backdiff.h"
#include "TempLat/lattice/algebra/spatialderivatives/forwdiff.h"
#include "TempLat/lattice/algebra/spatialderivatives/forwdij.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticeforwardgradient.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/algebra/spatialderivatives/neutdiff.h"
#include "TempLat/lattice/algebra/spatialderivatives/neutdij.h"
#include "TempLat/lattice/algebra/spatialderivatives/normgradientsquare.h"
#include "TempLat/lattice/algebra/listoperators/vectordotter.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/algebra/listoperators/listoperators.h"

namespace TempLat
{
  /** @brief The offsets at which expressions read fields, and the ghost depth check built on them. The stencils
   *  are read off onBeforeAssignment, which returns the smallest range containing 0 and every offset at which the
   *  right-hand side reads the target. */
  struct GhostReachTester {
    static void Test(TDDAssertion &tdd);
  };

  inline void GhostReachTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    const device::Idx nGrid = 8;
    auto toolBox1 = MemoryToolBox<NDim>::makeShared(nGrid, 1);
    auto toolBox2 = MemoryToolBox<NDim>::makeShared(nGrid, 2);

    const GhostReach none{};
    auto along = [](size_t axis, int lo, int hi) { return GhostReach{}.widened(axis, lo, hi); };
    auto all = [](int lo, int hi) { return GhostReach{}.widenedAll(NDim, lo, hi); };
    auto hull = [](GhostReach a, const GhostReach &b) { return a |= b; };

    // The range arithmetic itself.
    tdd.verify(none.isZero());
    tdd.verify(none.shifted(0, 1).shifted(1, 1).depth() == 1);
    tdd.verify(none.shifted(0, 1).shifted(0, -1).isZero());
    tdd.verify(along(2, -2, 1).depth() == 2);
    tdd.verify(hull(along(0, 0, 1), along(0, -2, 0)) == along(0, -2, 1));

    // Scalar stencils.
    {
      Field<double, NDim> phi("phi", toolBox2);
      Field<double, NDim> chi("chi", toolBox2);
      phi = 1.;
      chi = 1.;

      tdd.verify(phi.onBeforeAssignment(shift<1>(phi)) == along(0, 0, 1));
      tdd.verify(phi.onBeforeAssignment(shift<-3>(phi)) == along(2, -1, 0));
      tdd.verify(phi.onBeforeAssignment(shift<1, 0, -2>(phi)) == hull(along(0, 0, 1), along(2, -2, 0)));
      tdd.verify(phi.onBeforeAssignment(forwDiff(phi, 2_c)) == along(1, 0, 1));
      tdd.verify(phi.onBeforeAssignment(backDiff(phi, 2_c)) == along(1, -1, 0));
      tdd.verify(phi.onBeforeAssignment(neutDiff(phi, 2_c)) == along(1, -1, 1));
      tdd.verify(phi.onBeforeAssignment(forwDij(phi, 3_c)) == along(2, 0, 1));
      tdd.verify(phi.onBeforeAssignment(neutDij(phi, 3_c)) == along(2, -1, 1));
      tdd.verify(phi.onBeforeAssignment(LatLapl(phi)) == all(-1, 1));
      tdd.verify(phi.onBeforeAssignment(Grad2(phi)) == all(0, 1));
      tdd.verify(phi.onBeforeAssignment(dot(LatForwardGrad(phi), LatForwardGrad(chi))) == all(0, 1));

      // Nesting adds the stencils per axis; diagonal reads need one layer per axis.
      tdd.verify(phi.onBeforeAssignment(LatLapl(LatLapl(phi))) == all(-2, 2));
      tdd.verify(phi.onBeforeAssignment(shift<1>(shift<2>(phi))) == hull(along(0, 0, 1), along(1, 0, 1)));
      tdd.verify(phi.onBeforeAssignment(forwDiff(backDiff(phi, 1_c), 1_c)) == along(0, -1, 1));
      tdd.verify(phi.onBeforeAssignment(shift<1>(shift<-1>(phi))).isZero());

      // Branches: each read of the target counts, reads of other fields do not.
      tdd.verify(phi.onBeforeAssignment(phi + shift<2>(chi)).isZero());
      tdd.verify(phi.onBeforeAssignment(shift<1>(phi) + 2. * shift<-2>(phi)) == hull(along(0, 0, 1), along(1, -1, 0)));
      tdd.verify(phi.onBeforeAssignment(shift<1>(chi * LatLapl(phi))) == none.shifted(0, 1).widenedAll(NDim, -1, 1));
    }

    // Ghost depth check, nGhost = 1.
    {
      Field<double, NDim> phi("phi", toolBox1);
      Field<double, NDim> chi("chi", toolBox1);
      phi = 1.;
      chi = 1.;

      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = LatLapl(LatLapl(phi)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = shift<1>(shift<1>(phi)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = shift<0, 2, 0>(phi); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = forwDiff(forwDiff(phi, 1_c), 1_c); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = backDiff(backDiff(phi, 2_c), 2_c); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = neutDiff(neutDiff(phi, 3_c), 3_c); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = Grad2(Grad2(phi)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = phi + shift<1>(LatLapl(chi)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { average(LatLapl(LatLapl(phi))); }));

      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { chi = LatLapl(phi); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { chi = shift<1>(shift<2>(shift<3>(phi))); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { chi = forwDiff(backDiff(phi, 1_c), 1_c); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { chi = LatLapl(shift<1>(shift<-1>(phi))); }));

      // The message names the field.
      try {
        chi = LatLapl(LatLapl(phi));
      } catch (const GhostDepthTooSmall &e) {
        tdd.verify(std::string(e.what()).find("phi") != std::string::npos);
      }

      // A point-wise read needs no exchange, even through two opposite shifts.
      phi = 2. * phi;
      chi = shift<1>(shift<-1>(phi));
      tdd.verify(phi.getMemoryManager()->areGhostsStale());
      chi = shift<1>(phi);
      tdd.verify(!phi.getMemoryManager()->areGhostsStale());
    }

    // Composite fields, nGhost = 1.
    {
      SU2Field<double, NDim> U("U", toolBox1);
      SU2Field<double, NDim> V("V", toolBox1);
      ComplexField<double, NDim> psi("psi", toolBox1);
      SU2Doublet<double, NDim> D("D", toolBox1);
      SymTracelessField<double, NDim> h("h", toolBox1);

      tdd.verify(U.onBeforeAssignment(V * shift<1>(shift<2>(U))) == hull(along(0, 0, 1), along(1, 0, 1)));
      tdd.verify(U.onBeforeAssignment(shift<0, -1, 0>(U)) == along(1, -1, 0));
      tdd.verify(psi.onBeforeAssignment(shift<-2>(psi)) == along(1, -1, 0));
      tdd.verify(D.onBeforeAssignment(U * shift<3>(D)) == along(2, 0, 1));
      tdd.verify(h.onBeforeAssignment(shift<0, 0, 1>(h)) == along(2, 0, 1));

      tdd.verify(Throws<GhostDepthTooSmall>([&]() { V = shift<1>(shift<1>(U)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { D = U * shift<2>(shift<2>(D)); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { V = U * shift<1>(shift<2>(U)); }));

      ComplexField<double, NDim> chi("chi", toolBox1);
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = psi + shift<-3>(shift<-3>(psi)); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { chi = shift<3>(shift<-3>(shift<-3>(psi))); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { chi = psi * shift<0, 2, 0>(psi); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { chi = psi * shift<1>(shift<2>(psi)); }));

      SymTracelessField<double, NDim> g("g", toolBox1);
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { g = h + shift<-1>(shift<-1>(h)); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { g = h + shift<-1>(shift<2>(h)); }));

      SU2LieAlgebraField<double, NDim> A("A", toolBox1, LatticeParameters<double>());
      SU2LieAlgebraField<double, NDim> B("B", toolBox1, LatticeParameters<double>());
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { B = shift<2>(shift<2>(A)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { B = A + shift<0, 0, -2>(A); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { B = A + shift<1>(shift<-2>(A)); }));
    }

    // Matrix wrappers, which enter assignments through traces, nGhost = 1.
    {
      Field<double, NDim> phi("phi", toolBox1);
      Field<double, NDim> chi("chi", toolBox1);
      SymTracelessField<double, NDim> h("h", toolBox1);
      const auto far = shift<1>(shift<1>(phi));
      const auto near = shift<1>(phi);

      tdd.verify(
          Throws<GhostDepthTooSmall>([&]() { chi = multiplyTrace(ConstructSym(phi, phi, far, phi, phi, phi), h); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>(
          [&]() { chi = multiplyTrace(ConstructSym(phi, phi, near, phi, phi, phi), h); }));
      tdd.verify(Throws<GhostDepthTooSmall>(
          [&]() { chi = multiplyTrace(ConstructSym(phi, phi, phi, phi, phi, phi), shift<2>(shift<2>(h))); }));

      const auto c = Complexify(phi, phi);
      tdd.verify(
          Throws<GhostDepthTooSmall>([&]() { chi = Real(multiplyTrace(ConstructHerm(phi, c, c, far, c, phi), h)); }));
      tdd.verify(Throws<GhostDepthTooSmall>(
          [&]() { chi = Real(multiplyTrace(ConstructHerm(phi, c, Complexify(phi, far), phi, c, phi), h)); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>(
          [&]() { chi = Real(multiplyTrace(ConstructHerm(phi, c, c, near, c, phi), h)); }));

      const auto M = ConstructMatrix3x3(phi, phi, phi, phi, phi, phi, phi, phi, phi);
      tdd.verify(Throws<GhostDepthTooSmall>(
          [&]() { chi = multiplyTrace(M, ConstructMatrix3x3(phi, phi, phi, phi, phi, phi, phi, far, phi)); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>(
          [&]() { chi = multiplyTrace(M, ConstructMatrix3x3(phi, phi, phi, phi, phi, phi, phi, near, phi)); }));
    }

    // Field collections, nGhost = 1. They assign component by component.
    {
      FieldCollection<Field<double, NDim>, 2> fc("fc", toolBox1);
      FieldCollection<Field<double, NDim>, 2> gc("gc", toolBox1);

      tdd.verify(Throws<GhostDepthTooSmall>([&]() { gc = LatLapl(LatLapl(fc)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { gc = shift<2>(shift<2>(fc)); }));
      tdd.verify(Throws<GhostDepthTooSmall>([&]() { gc = fc + shift<0, -2, 0>(fc); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { gc = LatLapl(fc); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { gc = shift<1>(shift<2>(fc)); }));
      tdd.verify(DoesNotThrow<GhostDepthTooSmall>([&]() { gc = shift<1, -1, 1>(fc); }));
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GhostReachTester> test;
}
