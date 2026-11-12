/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/field/views/fieldviewconfig.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"

namespace TempLat
{
  /** @brief Detection of assignments whose right-hand side reads the target at a neighbouring site,
   *  as reported by onBeforeAssignment.
   *
   *  For the composite fields, the right-hand side is mostly a shifted copy of the target with permuted
   *  components (built from the same component fields). It reads every target component shifted, but never
   *  in the slot that component is assigned to.
   **/
  struct AssignmentAliasingTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace AssignmentAliasingDetail
  {
    /** @brief A scalar leaf counting how often the ghost walk reaches it. */
    struct CountingLeaf {
      int *calls;
      static consteval size_t getNDim() { return 3; }
      template <typename... IDX> DEVICE_INLINE_FUNCTION double eval(const IDX &...) const { return 0.; }
      device::Idx confirmGhostsUpToDate(const GhostReach &) const
      {
        ++*calls;
        return 0;
      }
    };
  } // namespace AssignmentAliasingDetail

  inline void AssignmentAliasingTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    const device::Idx nGrid = 8, nGhost = 2;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // Scalar fields.
    {
      Field<double, NDim> phi("phi", toolBox);
      Field<double, NDim> chi("chi", toolBox);
      phi = 1.;
      chi = 1.;

      // Point-wise self-reference is safe in place.
      tdd.verify(phi.onBeforeAssignment(phi).isZero());
      tdd.verify(phi.onBeforeAssignment(2. * phi + chi).isZero());

      // Non-local reads of other fields only.
      tdd.verify(phi.onBeforeAssignment(LatLapl(chi)).isZero());
      tdd.verify(phi.onBeforeAssignment(phi + shift<1>(chi)).isZero());

      // Non-local self-reference, directly and nested.
      tdd.verify(!phi.onBeforeAssignment(shift<1>(phi)).isZero());
      tdd.verify(!phi.onBeforeAssignment(2. * LatLapl(phi)).isZero());
      tdd.verify(!phi.onBeforeAssignment(phi + 0.1 * LatLapl(phi)).isZero());
      tdd.verify(!phi.onBeforeAssignment(LatLapl(LatLapl(phi))).isZero());
      tdd.verify(!phi.onBeforeAssignment(LatLapl(2. * phi + chi)).isZero());

      // A copy of the target shares its memory.
      Field<double, NDim> phiCopy = phi;
      tdd.verify(!phi.onBeforeAssignment(LatLapl(phiCopy)).isZero());

      // Detection does not depend on whether the ghosts actually need an exchange.
      phi = 2. * phi;
      tdd.verify(phi.getMemoryManager()->areGhostsStale());
      tdd.verify(!phi.onBeforeAssignment(LatLapl(phi)).isZero());

      // The target flag is cleared afterwards: phi as a plain source is not flagged.
      tdd.verify(chi.onBeforeAssignment(LatLapl(phi)).isZero());
    }

    // Complex fields.
    {
      ComplexField<double, NDim> psi("psi", toolBox);
      ComplexField<double, NDim> chi("chi", toolBox);
      ComplexField<double, NDim> psiSwapped(psi(1_c), psi(0_c));

      tdd.verify(psi.onBeforeAssignment(psi * psi).isZero());
      tdd.verify(psi.onBeforeAssignment(psi * shift<1>(chi)).isZero());
      tdd.verify(!psi.onBeforeAssignment(shift<1>(psi)).isZero());
      tdd.verify(!psi.onBeforeAssignment(shift<1>(psiSwapped)).isZero());
      tdd.verify(!psi.onBeforeAssignment(chi + 0.5 * shift<1>(psiSwapped)).isZero());
      tdd.verify(!psi.onBeforeAssignment(Complexify(shift<1>(Imag(psi)), Real(psi))).isZero());
    }

    // A composite assignment walks the right-hand side once, not once per target component.
    {
      ComplexField<double, NDim> psi("psi", toolBox);
      int calls = 0;
      AssignmentAliasingDetail::CountingLeaf leaf{&calls};
      psi.onBeforeAssignment(Complexify(leaf, leaf));
      tdd.verify(calls == 2);

      SU2Field<double, NDim> U("U", toolBox);
      calls = 0;
      U.onBeforeAssignment(leaf * U);
      tdd.verify(calls == 1);
    }

    // SU(2) fields.
    {
      SU2Field<double, NDim> U("U", toolBox);
      SU2Field<double, NDim> V("V", toolBox);

      tdd.verify(U.onBeforeAssignment(U * U).isZero());
      tdd.verify(U.onBeforeAssignment(shift<1>(V)).isZero());
      tdd.verify(!U.onBeforeAssignment(shift<1>(U)).isZero());
      tdd.verify(!U.onBeforeAssignment(V * dagger(shift<2>(U))).isZero());
    }

    // SU(2) doublets.
    {
      SU2Doublet<double, NDim> D("D", toolBox);
      SU2Doublet<double, NDim> E("E", toolBox);
      SU2Doublet<double, NDim> DPermuted(D(1_c), D(2_c), D(3_c), D(0_c));
      SU2Field<double, NDim> U("U", toolBox);
      ComplexField<double, NDim> psi("psi", toolBox);

      tdd.verify(D.onBeforeAssignment(U * D).isZero());
      tdd.verify(D.onBeforeAssignment(shift<1>(E)).isZero());
      tdd.verify(!D.onBeforeAssignment(shift<1>(D)).isZero());
      tdd.verify(!D.onBeforeAssignment(shift<1>(DPermuted)).isZero());
      tdd.verify(!D.onBeforeAssignment(U * shift<1>(DPermuted)).isZero());
      tdd.verify(!D.onBeforeAssignment(psi * shift<1>(DPermuted)).isZero());
    }

    // Symmetric traceless fields.
    {
      SymTracelessField<double, NDim> h("h", toolBox);
      SymTracelessField<double, NDim> g("g", toolBox);
      SymTracelessField<double, NDim> hPermuted(h(1_c), h(2_c), h(3_c), h(4_c), h(0_c));

      tdd.verify(h.onBeforeAssignment(2. * h).isZero());
      tdd.verify(h.onBeforeAssignment(h + shift<1>(g)).isZero());
      tdd.verify(!h.onBeforeAssignment(shift<1>(h)).isZero());
      tdd.verify(!h.onBeforeAssignment(shift<1>(hPermuted)).isZero());
      tdd.verify(!h.onBeforeAssignment(g + 2. * shift<1>(hPermuted)).isZero());
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::AssignmentAliasingTester> test;
}
