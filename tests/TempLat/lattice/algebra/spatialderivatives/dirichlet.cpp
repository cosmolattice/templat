
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

// Dirichlet-BC stencil tests. We choose a "bump" field whose value is exactly zero on the
// physical-domain boundary (global x[0] == 0 and global x[0] == nGrid-1 set to zero), and
// nonzero in the interior. Under Dirichlet BC the ghost slab is zero; combined with the
// boundary cells being zero, the Laplacian stencil at near-boundary interior cells is
// indistinguishable from a fully-interior stencil — i.e., no spurious boundary artifacts
// from the BC machinery. This is the canonical "open boundary" check.

#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/util/ndloop.h"
#include "TempLat/util/rangeiteration/tag.h"

namespace TempLat
{

  template <size_t NDim> struct DirichletStencilTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> inline void DirichletStencilTester<NDim>::Test(TDDAssertion &tdd)
  {
    const device::Idx nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();

    BCSpec<NDim> spec = allPeriodic<NDim>();
    spec[0] = BCType::Dirichlet;

    // Bump in dim 0: f(global_x[0] = k) = k * (nGrid - 1 - k). Zero at k=0 and k=nGrid-1
    // (so the field is zero on the physical-domain boundary in dim 0). Nonzero in interior.
    // Constant in other dims (Periodic). Centered Laplacian of this in dim 0:
    //   f[k+1] - 2 f[k] + f[k-1] = -2  for any interior k (it's a quadratic).
    // For the cells adjacent to the boundary (k=1 and k=nGrid-2), the stencil reaches f[0]
    // and f[nGrid-1] respectively — both ZERO under Dirichlet (and they ARE zero in the
    // interior data too). So the BC produces the same answer as the interior stencil here.
    Field<double, NDim> f("f_bump", toolBox, LatticeParameters<double>(), spec);
    SpatialCoordinate<NDim> x(toolBox);
    f = x(1_c) * (static_cast<double>(nGrid - 1) - x(1_c));
    f.updateGhosts();

    Field<double, NDim> lf("lf", toolBox);
    lf = LatLapl(f);

    // Other dims contribute 0 (constant in those dims, Periodic). Dim 0 contribution is -2
    // for every interior cell. The boundary cells themselves (k=0, k=nGrid-1) have value 0;
    // their stencil reaches into the Dirichlet ghost (= 0) on one side and an interior cell
    // on the other. At k=0:        f[-1]=0, f[0]=0, f[1] = 1*(nGrid-2) → laplacian = nGrid-2.
    // At k=nGrid-1: f[nGrid-1]=0, f[nGrid-2] = (nGrid-2)*1, ghost=0 → laplacian = nGrid-2.
    auto layout = toolBox->mLayouts.getConfigSpaceLayout();
    const auto &localSizes = layout.getLocalSizes();
    const auto &localStarts = layout.getLocalStarts();

    auto view = lf.getLocalNDHostView();
    bool ok = true;
    NDLoop<NDim>(view, [&](const auto &...indices) {
      device::IdxArray<NDim> idx{indices...};
      const ptrdiff_t globalK = localStarts[0] + idx[0];
      double expected = -2.0;
      if (globalK == 0 || globalK == nGrid - 1) expected = static_cast<double>(nGrid - 2);
      const double got = view(indices...);
      if (std::abs(got - expected) > 1e-12) {
        ok = false;
        std::stringstream ss;
        ss << "DirichletStencil mismatch at idx=" << idx << " globalK=" << globalK
           << " got=" << got << " expected=" << expected << " localSize0=" << localSizes[0] << "\n";
        sayMPI << ss.str();
      }
    });
    tdd.verify(ok);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::DirichletStencilTester<2>> bcTest2;
  TempLat::TDDContainer<TempLat::DirichletStencilTester<3>> bcTest3;
} // namespace
