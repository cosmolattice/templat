
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/parallel/device.h"
#include "TempLat/util/tdd/tdd.h"

#include "bctesthelpers.h"

namespace TempLat
{

  /** @brief Setter-path analogue of BCFillTester.
   *
   * Constructs a Field<double, NDim> with the 3-arg ctor (default all-periodic), then injects the
   * BC via setBCSpec() before assignment. Sweeps all 4 BC types per dimension and verifies the
   * low/high ghost slabs match what the ctor-path BCFillTester produces. Equivalence with the
   * ctor path is the contract: the setter is the same plumbing as the BC-aware ctor.
   */
  template <size_t NDim> struct FieldSetBCSpecTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void FieldSetBCSpecTester<NDim>::Test(TDDAssertion &tdd)
  {
    constexpr ptrdiff_t nGrid = 8;
    constexpr ptrdiff_t nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    for (size_t bcDim = 0; bcDim < NDim; ++bcDim) {
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Antiperiodic, nGrid, nGhost));
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Dirichlet, nGrid, nGhost));
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Neumann, nGrid, nGhost));
      tdd.verify(BCTestDetail::checkBCInDim<NDim>(toolBox, bcDim, BCType::Periodic, nGrid, nGhost));
    }
  }

  /** @brief Locks in the ghost-staleness contract on setBCSpec.
   *
   * Construct a Field with the default (all-periodic) BC. Assign a coordinate-dependent value and
   * call updateGhosts() — at this point ghosts are filled under the periodic BC. Then change the
   * BC to antiperiodic via setBCSpec() and call updateGhosts() again: the high-boundary ghost
   * row must now be sign-flipped. If setBCSpec failed to mark ghosts stale, the second
   * updateGhosts() would be a no-op and the ghost would still hold the periodic value.
   */
  template <size_t NDim> struct FieldSetBCSpecAfterAssignmentTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void FieldSetBCSpecAfterAssignmentTester<NDim>::Test(TDDAssertion &tdd)
  {
    constexpr ptrdiff_t nGrid = 8;
    constexpr ptrdiff_t nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    for (size_t bcDim = 0; bcDim < NDim; ++bcDim) {
      Field<double, NDim> f("f_setbc_after_assign", toolBox, LatticeParameters<double>());
      BCTestDetail::assignCoordinatePlusOne<NDim>(f, toolBox, bcDim);
      f.updateGhosts();

      BCSpec<NDim> spec = allPeriodic<NDim>();
      spec[bcDim] = BCType::Antiperiodic;
      f.setBCSpec(spec);
      f.updateGhosts();

      tdd.verify(BCTestDetail::verifyGhostFaces<NDim>(f, bcDim, BCType::Antiperiodic, nGrid, nGhost));
    }
  }

  /** @brief Locks in that BC mutation is genuinely live, not cached.
   *
   * Iterate through all four BC types in sequence on a single Field, calling setBCSpec() +
   * updateGhosts() between each, and verifying the ghost cells match the latest BC. A naive
   * implementation that cached the BC at first ghost-fill would fail every check after the
   * first.
   */
  template <size_t NDim> struct FieldSetBCSpecRoundTripTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void FieldSetBCSpecRoundTripTester<NDim>::Test(TDDAssertion &tdd)
  {
    constexpr ptrdiff_t nGrid = 8;
    constexpr ptrdiff_t nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    for (size_t bcDim = 0; bcDim < NDim; ++bcDim) {
      Field<double, NDim> f("f_setbc_roundtrip", toolBox, LatticeParameters<double>());
      BCTestDetail::assignCoordinatePlusOne<NDim>(f, toolBox, bcDim);

      for (BCType bc : {BCType::Antiperiodic, BCType::Dirichlet, BCType::Neumann, BCType::Periodic,
                        BCType::Antiperiodic}) {
        BCSpec<NDim> spec = allPeriodic<NDim>();
        spec[bcDim] = bc;
        f.setBCSpec(spec);
        // setBCSpec marks ghosts stale; the last assignment also did so. Either way, the next
        // updateGhosts() must refill under the new BC.
        f.updateGhosts();
        tdd.verify(BCTestDetail::verifyGhostFaces<NDim>(f, bcDim, bc, nGrid, nGhost));
      }
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::FieldSetBCSpecTester<2>> setBCSpecTest2;
  TempLat::TDDContainer<TempLat::FieldSetBCSpecTester<3>> setBCSpecTest3;

  TempLat::TDDContainer<TempLat::FieldSetBCSpecAfterAssignmentTester<2>> setBCSpecAfterAssignmentTest2;
  TempLat::TDDContainer<TempLat::FieldSetBCSpecAfterAssignmentTester<3>> setBCSpecAfterAssignmentTest3;

  TempLat::TDDContainer<TempLat::FieldSetBCSpecRoundTripTester<2>> setBCSpecRoundTripTest2;
  TempLat::TDDContainer<TempLat::FieldSetBCSpecRoundTripTester<3>> setBCSpecRoundTripTest3;
} // namespace
