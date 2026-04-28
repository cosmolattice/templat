#ifndef TEMPLAT_TESTS_GHOSTCELLS_BCTESTHELPERS_H
#define TEMPLAT_TESTS_GHOSTCELLS_BCTESTHELPERS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/util/ndloop.h"

#include <sstream>

namespace TempLat
{
  namespace BCTestDetail
  {
    inline double expectedLowGhost(BCType bc, ptrdiff_t globalSize, ptrdiff_t localStart, bool isLowBoundary)
    {
      if (!isLowBoundary) return static_cast<double>(localStart);
      switch (bc) {
      case BCType::Periodic:     return static_cast<double>(globalSize);
      case BCType::Antiperiodic: return -static_cast<double>(globalSize);
      case BCType::Dirichlet:    return 0.0;
      case BCType::Neumann:      return 1.0;
      }
      return 0.0;
    }

    inline double expectedHighGhost(BCType bc, ptrdiff_t globalSize, ptrdiff_t localStart, ptrdiff_t localSize,
                                    bool isHighBoundary)
    {
      if (!isHighBoundary) return static_cast<double>(localStart + localSize + 1);
      switch (bc) {
      case BCType::Periodic:     return 1.0;
      case BCType::Antiperiodic: return -1.0;
      case BCType::Dirichlet:    return 0.0;
      case BCType::Neumann:      return static_cast<double>(globalSize);
      }
      return 0.0;
    }

    enum class BCApplyMode { Constructor, Setter };

    template <size_t NDim>
    inline void assignCoordinatePlusOne(Field<double, NDim> &f,
                                        device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t bcDim)
    {
      SpatialCoordinate<NDim> x(toolBox);
      // Initialize so each cell holds (global_x[bcDim] + 1). SpatialCoordinate's selector is a
      // compile-time tag, so we dispatch on bcDim via constexpr_for.
      constexpr_for<1, NDim + 1>([&](auto dirTag) {
        constexpr size_t d = static_cast<size_t>(decltype(dirTag)::value) - 1;
        if (d == bcDim) f = x(dirTag) + 1.0;
      });
    }

    template <size_t NDim>
    bool verifyGhostFaces(const Field<double, NDim> &f, size_t bcDim, BCType bc, ptrdiff_t nGrid, ptrdiff_t nGhost)
    {
      auto toolBox = f.getToolBox();
      auto layout = toolBox->mLayouts.getConfigSpaceLayout();
      const auto &localSizes = layout.getLocalSizes();
      const auto &localStarts = layout.getLocalStarts();

      const bool isLowBoundary  = (localStarts[bcDim] == 0);
      const bool isHighBoundary = (localStarts[bcDim] + localSizes[bcDim] == nGrid);

      auto view = f.getFullNDHostView();

      auto walkFace = [&](bool low, bool isBoundary) {
        bool ok = true;
        const ptrdiff_t bcIdx = low ? (nGhost - 1) : (nGhost + localSizes[bcDim]);
        const double expected = low
            ? expectedLowGhost(bc, nGrid, localStarts[bcDim], isBoundary)
            : expectedHighGhost(bc, nGrid, localStarts[bcDim], localSizes[bcDim], isBoundary);

        device::IdxArray<NDim> idx{};
        for (size_t i = 0; i < NDim; ++i) idx[i] = (i == bcDim) ? bcIdx : nGhost;
        const ptrdiff_t loopDim = (bcDim == 0) ? 1 : 0;
        const ptrdiff_t loopExtent = (NDim == 1) ? 1 : localSizes[loopDim];
        for (ptrdiff_t i = 0; i < loopExtent; ++i) {
          if constexpr (NDim > 1) idx[loopDim] = nGhost + i;
          const double got = device::apply([&](auto... a) { return view(a...); }, idx);
          if (std::abs(got - expected) > 1e-14) {
            ok = false;
            std::stringstream ss;
            ss << "BCFill mismatch (bcDim=" << bcDim << ", bc=" << static_cast<int>(bc)
               << ", " << (low ? "low" : "high") << " face) at view_idx=" << idx
               << " got=" << got << " expected=" << expected
               << " localStart=" << localStarts[bcDim] << " localSize=" << localSizes[bcDim] << "\n";
            sayMPI << ss.str();
          }
        }
        return ok;
      };

      bool ok = true;
      ok &= walkFace(true,  isLowBoundary);
      ok &= walkFace(false, isHighBoundary);
      return ok;
    }

    /** @brief Construct a Field with the given BC (via ctor or via setter), set every cell to
     * (global_x[bcDim] + 1), updateGhosts(), and verify the low/high ghost slabs along bcDim.
     */
    template <size_t NDim>
    bool checkBCInDim(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t bcDim, BCType bc,
                      ptrdiff_t nGrid, ptrdiff_t nGhost, BCApplyMode mode = BCApplyMode::Constructor)
    {
      BCSpec<NDim> spec = allPeriodic<NDim>();
      spec[bcDim] = bc;

      auto runCheck = [&](Field<double, NDim> &f) {
        assignCoordinatePlusOne<NDim>(f, toolBox, bcDim);
        f.updateGhosts();
        return verifyGhostFaces<NDim>(f, bcDim, bc, nGrid, nGhost);
      };

      if (mode == BCApplyMode::Constructor) {
        Field<double, NDim> f("f_bc_check", toolBox, LatticeParameters<double>(), spec);
        return runCheck(f);
      } else {
        Field<double, NDim> f("f_bc_check", toolBox, LatticeParameters<double>());
        f.setBCSpec(spec);
        return runCheck(f);
      }
    }
  } // namespace BCTestDetail
} // namespace TempLat

#endif
