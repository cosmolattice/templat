
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

// Regression guard: the MPI Cartesian comm must remain all-periodic regardless of the user's
// BCSpec. The bd-551 design deliberately keeps Cart fully periodic (mPeriodic = all-1) so that
// neighbor-rank topology always wraps; non-Periodic BC is then applied as a post-step on
// boundary ranks (see GhostUpdater::applyLocalBCAtDimDepth). If someone later "fixes" this by
// plumbing BCSpec into MPICartesianGroup, the boundary-rank post-step logic in GhostUpdater
// would silently break (no exchange happens for non-periodic Cart axes). This test catches
// that regression by asking MPI itself for the comm's periodicity via MPI_Cart_get.

#include "TempLat/parallel/mpi/cartesian/mpicartesiangroup.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{
  struct CartesianPeriodicityTester {
    static void Test(TDDAssertion &tdd);
  };

  void CartesianPeriodicityTester::Test(TDDAssertion &tdd)
  {
#ifdef HAVE_MPI
    constexpr size_t NDim = 3;
    constexpr device::Idx nGrid = 8, nGhost = 1;

    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    // Ask MPI for the comm's periodicity. Even if internal accessors gain a setter or default
    // changes, MPI_Cart_get returns ground truth from the comm itself.
    int dims[NDim], periods[NDim], coords[NDim];
    int rc = MPI_Cart_get(toolBox->mGroup.getComm(), static_cast<int>(NDim), dims, periods, coords);
    tdd.verify(rc == MPI_SUCCESS);

    bool allPeriodic = true;
    for (size_t i = 0; i < NDim; ++i) allPeriodic = allPeriodic && (periods[i] == 1);
    tdd.verify(allPeriodic);
#else
    // No-op under non-MPI builds — there's no Cart comm to inspect. Single tdd assertion keeps
    // the test discoverable in summaries.
    tdd.verify(true);
#endif
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::CartesianPeriodicityTester> test;
} // namespace
