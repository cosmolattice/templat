
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026
//
// Per-component boundary conditions through the COALESCED ghost path.
//
// Every multi-component field (ComplexField, SU2Field, SU2Doublet, SU2LieAlgebraField,
// SymTracelessField) updates its ghosts through MemoryManager::updateGhostsBatch /
// confirmGhostsUpToDateBatch, which fuse the C component blocks into one exchange. That batch
// originally carried ONE BCSpec for the whole field and rejected any batch whose components
// disagreed:
//
//     "confirmGhostsUpToDateBatch: components do not share a BCSpec. All components of a
//      multi-component field must carry the same boundary condition."
//
// That invariant is false. C-star boundary conditions are exactly a per-component BC: the
// identification Phi(x + L) = Phi*(x) is, in the real-component basis these fields are stored
// in, a sign flip on some components and not others --
//     (c0, c1, c2, c3) -> (c0, -c1, c2, -c3)   for an SU(2) group element and for an SU(2) doublet,
//     (Re, Im)         -> (Re, -Im)            for a complex scalar.
// So a legitimate batch mixes Periodic and Antiperiodic, and the ghost slab of an antiperiodic
// component must hold MINUS its global-wrap partner while a periodic sibling holds PLUS it, in
// the same dimension, in the same launch.
//
// The tests below drive the batch entry points with mixed BCSpecs and assert the resulting ghost
// contents per component, per dimension. Against a build that still carries the uniformity
// invariant every one of them fails on the throw above.
//
// Companion: ghostupdaterpercomponentbc.cpp covers GhostUpdater::updateBatch's per-component
// span overload directly, one level below MemoryManager.

#include "TempLat/lattice/algebra/complexalgebra/complexfield.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfield.h"
#include "TempLat/lattice/algebra/su2algebra/su2doublet.h"
#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/algebra/su2algebra/su2liealgebrafield.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/memory/memorymanager.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/parallel/device.h"
#include "TempLat/util/tdd/tdd.h"

#include "bctesthelpers.h"

#include <array>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace TempLat
{
  namespace PerComponentBCDetail
  {
    /** @brief Runs a check that is expected to complete, turning any escaping exception into a
     *  reported failure. Pre-fix builds throw out of the batch entry point; without this the whole
     *  test binary would die on the first assertion and report nothing useful. */
    template <typename F> bool guarded(const char *what, F &&fn)
    {
      try {
        return fn();
      } catch (const std::exception &e) {
        std::stringstream ss;
        ss << "[" << what << "] batch path threw instead of honouring the per-component BCSpecs: " << e.what() << "\n";
        sayMPI << ss.str();
        return false;
      }
    }

    template <size_t NDim> BCSpec<NDim> specWith(size_t bcDim, BCType bc)
    {
      BCSpec<NDim> s = allPeriodic<NDim>();
      s[bcDim] = bc;
      return s;
    }

    /** @brief The C-star pattern of a 4-real-component SU(2) object: c1 and c3 flip sign under
     *  complex conjugation, c0 and c2 do not. */
    inline std::vector<BCType> cStarSU2(BCType flip = BCType::Antiperiodic)
    {
      return {BCType::Periodic, flip, BCType::Periodic, flip};
    }

    /** @brief Fill with (global_x[bcDim] + 1) so every BC has a distinct, analytically known ghost
     *  value -- see bctesthelpers.h. */
    template <size_t NDim>
    void fillCoord(Field<double, NDim> &f, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t bcDim)
    {
      BCTestDetail::assignCoordinatePlusOne<NDim>(f, toolBox, bcDim);
    }

    /** @brief Fill with a value that depends on EVERY dimension and on the component index. A
     *  fill that varies in one dimension only cannot see a corner ghost taken from the wrong
     *  dimension, and identical values across components cannot see a coalesced pack that swaps
     *  two components' slabs. */
    template <size_t NDim>
    void fillAllDimsLabelled(Field<double, NDim> &f, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t c)
    {
      SpatialCoordinate<NDim> x(toolBox);
      const double label = 1.0e6 * static_cast<double>(c + 1);
      static_assert(NDim == 2 || NDim == 3, "fillAllDimsLabelled covers the dimensions the tests use.");
      if constexpr (NDim == 2)
        f = x(Tag<1>()) + x(Tag<2>()) * 37.0 + (label + 1.0);
      else if constexpr (NDim == 3)
        f = x(Tag<1>()) + x(Tag<2>()) * 37.0 + x(Tag<3>()) * 1369.0 + (label + 1.0);
    }

    template <size_t NDim>
    size_t paddedCellCount(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t nGhost)
    {
      const auto &localSizes = toolBox->mLayouts.getConfigSpaceLayout().getLocalSizes();
      size_t total = 1;
      for (size_t i = 0; i < NDim; ++i)
        total *= static_cast<size_t>(localSizes[i]) + 2 * nGhost;
      return total;
    }

    /** @brief Drive C scalar Fields carrying per-component BCSpecs through the batch entry point and
     *  assert each component's ghost slabs hold what ITS OWN boundary condition prescribes.
     *
     *  This is the load-bearing check: with `bcs = {Periodic, Antiperiodic, Periodic, Antiperiodic}`
     *  components 1 and 3 must come back sign-flipped at the global boundary while 0 and 2 must not,
     *  out of a single fused update of the same dimension. `useConfirm` picks
     *  confirmGhostsUpToDateBatch (the entry point that fired in production) over updateGhostsBatch. */
    template <size_t NDim>
    bool checkMixedBatch(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t bcDim,
                         const std::vector<BCType> &bcs, ptrdiff_t nGrid, ptrdiff_t nGhost, bool useConfirm)
    {
      const size_t C = bcs.size();
      std::vector<std::unique_ptr<Field<double, NDim>>> fields;
      std::vector<MemoryManager<double, NDim> *> mgrs;
      fields.reserve(C);
      mgrs.reserve(C);
      for (size_t c = 0; c < C; ++c) {
        fields.push_back(std::make_unique<Field<double, NDim>>("pcbc_" + std::to_string(c), toolBox,
                                                               LatticeParameters<double>()));
        fields[c]->setBCSpec(specWith<NDim>(bcDim, bcs[c]));
        fillCoord<NDim>(*fields[c], toolBox, bcDim);
        mgrs.push_back(fields[c]->getMemoryManager().get());
      }

      auto span = std::span<MemoryManager<double, NDim> *const>(mgrs.data(), mgrs.size());
      if (useConfirm)
        MemoryManager<double, NDim>::confirmGhostsUpToDateBatch(span);
      else
        MemoryManager<double, NDim>::updateGhostsBatch(span);

      bool ok = true;
      for (size_t c = 0; c < C; ++c)
        ok &= BCTestDetail::verifyGhostFaces<NDim>(*fields[c], bcDim, bcs[c], nGrid, nGhost);

      // confirmGhostsUpToDateBatch only exchanges the components whose ghosts are STALE, and it has
      // to keep its per-component BCSpec list aligned with that filtered block list. Dirty a single
      // antiperiodic component and re-confirm: a BC list indexed against the unfiltered component
      // list would hand component 1's blocks component 0's (periodic) BC here.
      if (useConfirm && C > 1) {
        const size_t dirty = 1;
        fillCoord<NDim>(*fields[dirty], toolBox, bcDim);
        MemoryManager<double, NDim>::confirmGhostsUpToDateBatch(span);
        ok &= BCTestDetail::verifyGhostFaces<NDim>(*fields[dirty], bcDim, bcs[dirty], nGrid, nGhost);
      }
      return ok;
    }

    /** @brief The coalesced batch must be indistinguishable, cell for cell, from updating each
     *  component on its own with that component's BCSpec -- ghosts, corners and edges included.
     *
     *  Per-component `Field::updateGhosts()` takes the single-block GhostUpdater path, which is the
     *  path all the other BC tests in this tree cover; the batch takes the fused/coalesced one. Any
     *  divergence -- a mis-partitioned fused kernel launch, a BC applied to the wrong component, a
     *  corner filled from the wrong dimension -- shows up as a mismatch. Comparison is exact: this
     *  is data movement and negation, not arithmetic, so "agrees" means "bit for bit". */
    template <size_t NDim>
    bool checkBatchMatchesSingle(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
                                 const std::vector<BCSpec<NDim>> &specs, size_t nGhost)
    {
      const size_t C = specs.size();
      const size_t total = paddedCellCount<NDim>(toolBox, nGhost);

      std::vector<std::unique_ptr<Field<double, NDim>>> batched, single;
      std::vector<MemoryManager<double, NDim> *> mgrs;
      batched.reserve(C);
      single.reserve(C);
      mgrs.reserve(C);
      for (size_t c = 0; c < C; ++c) {
        batched.push_back(std::make_unique<Field<double, NDim>>("bvs_b_" + std::to_string(c), toolBox,
                                                                LatticeParameters<double>()));
        single.push_back(std::make_unique<Field<double, NDim>>("bvs_s_" + std::to_string(c), toolBox,
                                                               LatticeParameters<double>()));
        batched[c]->setBCSpec(specs[c]);
        single[c]->setBCSpec(specs[c]);
        fillAllDimsLabelled<NDim>(*batched[c], toolBox, c);
        fillAllDimsLabelled<NDim>(*single[c], toolBox, c);
        mgrs.push_back(batched[c]->getMemoryManager().get());
      }

      MemoryManager<double, NDim>::updateGhostsBatch(
          std::span<MemoryManager<double, NDim> *const>(mgrs.data(), mgrs.size()));
      for (size_t c = 0; c < C; ++c)
        single[c]->updateGhosts();

      bool ok = true;
      for (size_t c = 0; c < C; ++c) {
        auto vb = batched[c]->getRawHostView();
        auto vs = single[c]->getRawHostView();
        for (size_t i = 0; i < total; ++i) {
          if (vb(i) != vs(i)) {
            ok = false;
            std::stringstream ss;
            ss << "batch/single mismatch: NDim=" << NDim << " nGhost=" << nGhost << " C=" << C << " comp=" << c
               << " flat=" << i << " batch=" << vb(i) << " single=" << vs(i) << "\n";
            sayMPI << ss.str();
            break; // one report per component is enough to localize it
          }
        }
      }
      return ok;
    }

    /** @brief Set a per-component BC on a real multi-component field type, fill, update its ghosts
     *  through whichever batch entry point the type uses, and verify every component's ghost slabs.
     *  `Cs...` are the component tags of the type (SU2LieAlgebraField's start at 1; its Tag<0> is a
     *  compile-time zero, not a Field). */
    template <size_t NDim, typename FieldT, int... Cs>
    bool checkMixedFieldType(FieldT &fld, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox, size_t bcDim,
                             const std::vector<BCType> &bcs, ptrdiff_t nGrid, ptrdiff_t nGhost, bool useConfirm,
                             std::integer_sequence<int, Cs...>)
    {
      static_assert(sizeof...(Cs) > 1, "A per-component BC batch needs more than one component.");
      size_t k = 0;
      (
          [&] {
            fld(Tag<Cs>()).setBCSpec(specWith<NDim>(bcDim, bcs[k]));
            fillCoord<NDim>(fld(Tag<Cs>()), toolBox, bcDim);
            ++k;
          }(),
          ...);

      if (useConfirm) {
        if constexpr (requires { fld.confirmGhostsUpToDate(); })
          fld.confirmGhostsUpToDate();
        else
          fld.updateGhosts();
      } else {
        fld.updateGhosts();
      }

      bool ok = true;
      k = 0;
      (
          [&] {
            ok &= BCTestDetail::verifyGhostFaces<NDim>(fld(Tag<Cs>()), bcDim, bcs[k], nGrid, nGhost);
            ++k;
          }(),
          ...);
      return ok;
    }
  } // namespace PerComponentBCDetail

  /** @brief MemoryManager's batch entry points must accept, and honour, a per-component BCSpec.
   *
   * Sweeps every dimension. Patterns: the SU(2)/doublet C-star pattern (P, AP, P, AP), the complex
   * scalar's (Re, Im) -> (Re, -Im), a 4-component batch carrying all four distinct BCTypes at once
   * (the fused local kernel partitions by BCType-in-dim and launches once per distinct type, so
   * this is its widest partition), and a batch wider than the fused kernel's capture budget, which
   * falls back to the per-component path and must still read each component's own BC.
   *
   * Both entry points are covered: updateGhostsBatch (unconditional) and confirmGhostsUpToDateBatch
   * (stale-filtered) -- the latter is the one that aborted in production, and its BCSpec list is
   * built alongside the filtered block list, so it also gets a partial-staleness check.
   */
  template <size_t NDim> struct MemoryManagerPerComponentBCTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void MemoryManagerPerComponentBCTester<NDim>::Test(TDDAssertion &tdd)
  {
    using namespace PerComponentBCDetail;
    constexpr ptrdiff_t nGrid = 8;

    const std::vector<std::pair<const char *, std::vector<BCType>>> patterns = {
        {"su2-cstar", cStarSU2()},
        {"complex-cstar", {BCType::Periodic, BCType::Antiperiodic}},
        {"all-four-bctypes", {BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet, BCType::Neumann}},
        {"wide-batch-fallback",
         {BCType::Periodic, BCType::Antiperiodic, BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet,
          BCType::Neumann, BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet}},
    };

    for (ptrdiff_t nGhost : {1, 2}) {
      auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
      toolBox->unsetVerbose();
      for (size_t bcDim = 0; bcDim < NDim; ++bcDim) {
        for (const auto &[name, bcs] : patterns) {
          tdd.verify(guarded(name, [&] { return checkMixedBatch<NDim>(toolBox, bcDim, bcs, nGrid, nGhost, false); }));
          tdd.verify(guarded(name, [&] { return checkMixedBatch<NDim>(toolBox, bcDim, bcs, nGrid, nGhost, true); }));
        }
      }
    }
  }

  /** @brief The batch result must equal the per-component single-block result cell for cell.
   *
   * The uniform cases come first and are the fast-path regression guard: the per-component plumbing
   * added an anyNonPeriodic() early-out and a uniformity check in front of the fused kernel, and a
   * uniform batch must still come out exactly as it did before -- one BCType, unchanged behaviour.
   * The mixed cases are the new contract, including specs that differ in DIFFERENT dimensions, since
   * the partition is per-dimension and each dimension is decided on its own.
   */
  template <size_t NDim> struct BatchVsSinglePerComponentBCTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void BatchVsSinglePerComponentBCTester<NDim>::Test(TDDAssertion &tdd)
  {
    using namespace PerComponentBCDetail;
    constexpr ptrdiff_t nGrid = 8;

    for (size_t nGhost : {1u, 2u}) {
      auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, static_cast<ptrdiff_t>(nGhost));
      toolBox->unsetVerbose();

      // Uniform batches: unchanged behaviour required.
      for (auto bc : {BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet, BCType::Neumann}) {
        for (size_t d = 0; d < NDim; ++d) {
          std::vector<BCSpec<NDim>> uniform(4, specWith<NDim>(d, bc));
          tdd.verify(guarded("uniform", [&] { return checkBatchMatchesSingle<NDim>(toolBox, uniform, nGhost); }));
        }
      }

      // Mixed within one dimension: the C-star pattern, and all four BCTypes at once.
      for (size_t d = 0; d < NDim; ++d) {
        std::vector<BCSpec<NDim>> cstar;
        for (auto bc : cStarSU2())
          cstar.push_back(specWith<NDim>(d, bc));
        tdd.verify(guarded("mixed-cstar", [&] { return checkBatchMatchesSingle<NDim>(toolBox, cstar, nGhost); }));

        std::vector<BCSpec<NDim>> allFour;
        for (auto bc : {BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet, BCType::Neumann})
          allFour.push_back(specWith<NDim>(d, bc));
        tdd.verify(guarded("mixed-allfour", [&] { return checkBatchMatchesSingle<NDim>(toolBox, allFour, nGhost); }));
      }

      // Mixed across dimensions: component 1 is antiperiodic in dim 0 only, component 2 Dirichlet in
      // dim 1 only, component 3 mixes both. Each dimension must partition on its own.
      {
        std::vector<BCSpec<NDim>> crossDim(4, allPeriodic<NDim>());
        crossDim[1][0] = BCType::Antiperiodic;
        crossDim[2][NDim - 1] = BCType::Dirichlet;
        crossDim[3][0] = BCType::Antiperiodic;
        crossDim[3][NDim - 1] = BCType::Neumann;
        tdd.verify(guarded("mixed-crossdim", [&] { return checkBatchMatchesSingle<NDim>(toolBox, crossDim, nGhost); }));
      }
    }
  }

  /** @brief The real multi-component field types under C-star.
   *
   * SU2Field and SU2Doublet: (c0, c1, c2, c3) -> (c0, -c1, c2, -c3). ComplexField: (Re, Im) ->
   * (Re, -Im). SU2LieAlgebraField: for anti-hermitian A = i a_k sigma_k, A* flips a1 and a3.
   * SymTracelessField carries five real components and batches them the same way; it has no C-star
   * meaning of its own, so it is driven with a mixed pattern purely to pin the batch plumbing.
   *
   * Each of these types was latent behind the same uniformity invariant -- these are the call sites
   * that would have aborted the moment their model set a per-component BC.
   */
  template <size_t NDim> struct MultiComponentFieldPerComponentBCTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void MultiComponentFieldPerComponentBCTester<NDim>::Test(TDDAssertion &tdd)
  {
    using namespace PerComponentBCDetail;
    constexpr ptrdiff_t nGrid = 8;
    constexpr ptrdiff_t nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
    toolBox->unsetVerbose();
    LatticeParameters<double> latPar;

    const std::vector<BCType> su2CStar = cStarSU2();
    const std::vector<BCType> complexCStar = {BCType::Periodic, BCType::Antiperiodic};
    const std::vector<BCType> algebraCStar = {BCType::Antiperiodic, BCType::Periodic, BCType::Antiperiodic};
    const std::vector<BCType> symMixed = {BCType::Periodic, BCType::Antiperiodic, BCType::Periodic,
                                          BCType::Antiperiodic, BCType::Periodic};

    for (size_t bcDim = 0; bcDim < NDim; ++bcDim) {
      // SU2Field: updateGhosts() AND confirmGhostsUpToDate() -- the latter is the exact call that
      // aborted (SU2Field::confirmGhostsUpToDate -> confirmGhostsUpToDateBatch -> throw).
      tdd.verify(guarded("SU2Field", [&] {
        SU2Field<double, NDim> u("pcbc_su2", toolBox, latPar);
        return checkMixedFieldType<NDim>(u, toolBox, bcDim, su2CStar, nGrid, nGhost, false,
                                         std::integer_sequence<int, 0, 1, 2, 3>{});
      }));
      tdd.verify(guarded("SU2Field-confirm", [&] {
        SU2Field<double, NDim> u("pcbc_su2c", toolBox, latPar);
        return checkMixedFieldType<NDim>(u, toolBox, bcDim, su2CStar, nGrid, nGhost, true,
                                         std::integer_sequence<int, 0, 1, 2, 3>{});
      }));

      tdd.verify(guarded("SU2Doublet", [&] {
        SU2Doublet<double, NDim> h("pcbc_doublet", toolBox, latPar);
        return checkMixedFieldType<NDim>(h, toolBox, bcDim, su2CStar, nGrid, nGhost, false,
                                         std::integer_sequence<int, 0, 1, 2, 3>{});
      }));

      tdd.verify(guarded("ComplexField", [&] {
        ComplexField<double, NDim> z("pcbc_complex", toolBox, latPar);
        return checkMixedFieldType<NDim>(z, toolBox, bcDim, complexCStar, nGrid, nGhost, false,
                                         std::integer_sequence<int, 0, 1>{});
      }));

      tdd.verify(guarded("SU2LieAlgebraField", [&] {
        SU2LieAlgebraField<double, NDim> e("pcbc_algebra", toolBox, latPar);
        return checkMixedFieldType<NDim>(e, toolBox, bcDim, algebraCStar, nGrid, nGhost, true,
                                         std::integer_sequence<int, 1, 2, 3>{});
      }));

      tdd.verify(guarded("SymTracelessField", [&] {
        SymTracelessField<double, NDim> s("pcbc_symtl", toolBox, latPar);
        return checkMixedFieldType<NDim>(s, toolBox, bcDim, symMixed, nGrid, nGhost, false,
                                         std::integer_sequence<int, 0, 1, 2, 3, 4>{});
      }));
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MemoryManagerPerComponentBCTester<2>> mmPerComponentBC2;
  TempLat::TDDContainer<TempLat::MemoryManagerPerComponentBCTester<3>> mmPerComponentBC3;

  TempLat::TDDContainer<TempLat::BatchVsSinglePerComponentBCTester<2>> batchVsSingle2;
  TempLat::TDDContainer<TempLat::BatchVsSinglePerComponentBCTester<3>> batchVsSingle3;

  TempLat::TDDContainer<TempLat::MultiComponentFieldPerComponentBCTester<2>> multiComponent2;
  TempLat::TDDContainer<TempLat::MultiComponentFieldPerComponentBCTester<3>> multiComponent3;
} // namespace
