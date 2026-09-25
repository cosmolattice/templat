#ifndef TEMPLAT_TESTS_GHOSTCELLS_BCTESTHELPERS_H
#define TEMPLAT_TESTS_GHOSTCELLS_BCTESTHELPERS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/util/ndloop.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

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

    /** @brief The sign a cell picks up from ONE direction.
     *
     * `crossed` says the ghost cell's source is the global-wrap partner in that direction, i.e.
     * the boundary condition applies there rather than an ordinary neighbour exchange.
     */
    inline double directionSign(BCType bc, bool crossed)
    {
      if (!crossed) return 1.0;
      return bc == BCType::Antiperiodic ? -1.0 : 1.0;
    }

    /** @brief The unique code every corner check fills the lattice with.
     *
     * `1 + x_0 + nGrid x_1 + nGrid^2 x_2`, i.e. the global position written in base nGrid. Two
     * properties matter and neither is negotiable:
     *
     *  - it is INJECTIVE over the global lattice, so a cell whose value is right can only have come
     *    from the one cell it should have come from. A fill that varies in a single dimension (what
     *    every other BC test in this tree uses) cannot tell a corner taken from the wrong dimension
     *    from a corner taken from the right one -- they hold the same number.
     *  - it is strictly positive, so a wrong SIGN is always visible, including at the origin.
     */
    template <size_t NDim> inline double globalCode(const std::array<ptrdiff_t, NDim> &g, ptrdiff_t nGrid)
    {
      double code = 1.0;
      double mult = 1.0;
      for (size_t d = 0; d < NDim; ++d) {
        code += static_cast<double>(g[d]) * mult;
        mult *= static_cast<double>(nGrid);
      }
      return code;
    }

    template <size_t NDim>
    inline void assignGlobalCode(Field<double, NDim> &f, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
                                 ptrdiff_t nGrid)
    {
      static_assert(NDim >= 1 && NDim <= 3, "assignGlobalCode covers the dimensions the BC tests use.");
      SpatialCoordinate<NDim> x(toolBox);
      const double N = static_cast<double>(nGrid);
      if constexpr (NDim == 1)
        f = x(Tag<1>()) + 1.0;
      else if constexpr (NDim == 2)
        f = x(Tag<1>()) + x(Tag<2>()) * N + 1.0;
      else
        f = x(Tag<1>()) + x(Tag<2>()) * N + x(Tag<3>()) * (N * N) + 1.0;
    }

    /** @brief Verify EVERY padded cell -- faces, edges and CORNERS -- against the product of the
     *  per-direction signs.
     *
     * Why this exists: verifyGhostFaces (above) pins one ghost index in ONE dimension and holds
     * every other index at `nGhost`, an interior cell. It therefore never reads a cell that is a
     * ghost in two directions at once. Before this function, NO test in either tree pinned a corner
     * value under any non-periodic boundary condition -- the whole edge/corner region was covered
     * only by "the batch agrees with the single-block path", which is satisfied just as well by two
     * identically wrong answers.
     *
     * The contract being pinned is the one a mixed spec makes non-trivial: the dimension sweep is
     * sequential, each dimension applies ITS OWN boundary condition to a slab that may already be a
     * ghost from an earlier dimension, so a cell that wrapped the global boundary in k directions
     * carries the product of those k directions' signs -- and NOT, say, the first one's, or the
     * last one's, or one sign for the whole cell.
     *
     * Only Periodic and Antiperiodic are handled: Dirichlet and Neumann overwrite rather than sign,
     * so "product of signs" is not their semantics and a spec containing one is rejected outright
     * rather than silently mis-verified.
     *
     * @param view    a full padded ND host view of the data.
     * @param layout  the config-space layout the view was built from.
     * @param spec    the per-direction BC the data was updated with.
     * @param value   value(globalCoords) -> the owned value at those GLOBAL coordinates.
     */
    template <size_t NDim, typename ViewT, typename LayoutT, typename ValueFn>
    bool verifyGhostCornersView(const ViewT &view, const LayoutT &layout, const BCSpec<NDim> &spec, ptrdiff_t nGrid,
                                ptrdiff_t nGhost, ValueFn &&value, const std::string &what)
    {
      for (size_t d = 0; d < NDim; ++d) {
        if (spec[d] != BCType::Periodic && spec[d] != BCType::Antiperiodic)
          throw std::logic_error("verifyGhostCorners: only Periodic/Antiperiodic have a product-of-signs "
                                 "semantics; Dirichlet and Neumann overwrite the ghost instead.");
      }

      const auto &localSizes = layout.getLocalSizes();
      const auto &localStarts = layout.getLocalStarts();

      device::IdxArray<NDim> extents{};
      size_t total = 1;
      for (size_t d = 0; d < NDim; ++d) {
        extents[d] = static_cast<device::Idx>(localSizes[d] + 2 * nGhost);
        total *= static_cast<size_t>(extents[d]);
      }

      bool ok = true;
      size_t reported = 0;
      device::IdxArray<NDim> idx{};
      for (size_t flat = 0; flat < total; ++flat) {
        // Odometer over the whole padded volume.
        size_t rest = flat;
        for (size_t d = 0; d < NDim; ++d) {
          idx[d] = static_cast<device::Idx>(rest % static_cast<size_t>(extents[d]));
          rest /= static_cast<size_t>(extents[d]);
        }

        std::array<ptrdiff_t, NDim> wrapped{};
        double sign = 1.0;
        for (size_t d = 0; d < NDim; ++d) {
          const ptrdiff_t g = localStarts[d] + static_cast<ptrdiff_t>(idx[d]) - nGhost;
          // "Crossed" means the GLOBAL boundary, not a rank boundary: on an interior rank a ghost
          // in a split dimension is an ordinary neighbour copy and carries no sign.
          const bool crossed = (g < 0) || (g >= nGrid);
          wrapped[d] = ((g % nGrid) + nGrid) % nGrid;
          sign *= directionSign(spec[d], crossed);
        }

        const double expected = sign * value(wrapped);
        const double got = device::apply([&](auto... a) { return view(a...); }, idx);
        // Data movement and negation only -- no arithmetic -- and every value here is an
        // exactly representable integer, so this is effectively an equality test. It has to be
        // tight in RELATIVE terms: the batch fill carries a ~1e9 component label, and a loose
        // relative tolerance there would swallow an off-by-one in the coordinate part, i.e.
        // exactly the "corner taken from the neighbouring cell" bug this function is for.
        if (std::abs(got - expected) > 1e-12 * std::max(1.0, std::abs(expected))) {
          ok = false;
          if (reported++ < 8) {
            std::stringstream ss;
            ss << "ghost-corner mismatch [" << what << "] at view_idx=" << idx << " got=" << got
               << " expected=" << expected << " (sign=" << sign << ")\n";
            sayMPI << ss.str();
          }
        }
      }
      return ok;
    }

    /** @brief verifyGhostCornersView for a Field filled by assignGlobalCode. */
    template <size_t NDim>
    bool verifyGhostCorners(const Field<double, NDim> &f, const BCSpec<NDim> &spec, ptrdiff_t nGrid, ptrdiff_t nGhost,
                            const std::string &what)
    {
      auto layout = f.getToolBox()->mLayouts.getConfigSpaceLayout();
      return verifyGhostCornersView<NDim>(
          f.getFullNDHostView(), layout, spec, nGrid, nGhost,
          [&](const std::array<ptrdiff_t, NDim> &g) { return globalCode<NDim>(g, nGrid); }, what);
    }

    /** @brief Build a Field with `spec`, fill it with the injective global code, updateGhosts(),
     *  and check the WHOLE padded volume including edges and corners.
     */
    template <size_t NDim>
    bool checkCornersForSpec(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, const BCSpec<NDim> &spec,
                             ptrdiff_t nGrid, ptrdiff_t nGhost, const std::string &what)
    {
      Field<double, NDim> f("f_corner_check", toolBox, LatticeParameters<double>());
      f.setBCSpec(spec);
      assignGlobalCode<NDim>(f, toolBox, nGrid);
      f.updateGhosts();
      return verifyGhostCorners<NDim>(f, spec, nGrid, nGhost, what);
    }

    /** @brief Construct a Field with the given BC, set every cell to
     * (global_x[bcDim] + 1), updateGhosts(), and verify the low/high ghost slabs along bcDim.
     */
    template <size_t NDim>
    bool checkBCInDim(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t bcDim, BCType bc,
                      ptrdiff_t nGrid, ptrdiff_t nGhost)
    {
      BCSpec<NDim> spec = allPeriodic<NDim>();
      spec[bcDim] = bc;
      Field<double, NDim> f("f_bc_check", toolBox, LatticeParameters<double>());
      f.setBCSpec(spec);
      assignCoordinatePlusOne<NDim>(f, toolBox, bcDim);
      f.updateGhosts();
      return verifyGhostFaces<NDim>(f, bcDim, bc, nGrid, nGhost);
    }
  } // namespace BCTestDetail
} // namespace TempLat

#endif
