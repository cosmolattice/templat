#ifndef TEMPLAT_LATTICE_MEASUREMENTS_RADIALPROJECTOR_H
#define TEMPLAT_LATTICE_MEASUREMENTS_RADIALPROJECTOR_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019

#include <algorithm>
#include <cstddef>

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionresult.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglequantity.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialbincomputer.h"
#include "TempLat/lattice/algebra/operators/squareroot.h"

#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojector_scratch.h"
#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/algebra/spacestateinterface.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"

namespace TempLat
{
  /** @brief A class which projects any N-D lattice on its positive-definite radial coordinate. In other
   *  words, integrating out all the angular dimensions.
   *  When in Fourier-space, this routine takes into account what the redundancies are.
   *  Only the unique values are counted. Real-valued entries get weight 0.5
   *  (they only contribute one float value out of the two for each complex value),
   *  and of the hermitian pairs, only the positivePartner is taken into account.
   *
   *  See `RadialProjectionResult` for a description of the result type.
   *
   * Unit test: ctest -R test-radialprojector
   **/

  template <typename T> class RadialProjector
  {
  public:
    using vType = typename GetGetReturnType<T>::type;
    using sType = typename GetFloatType<vType>::type;

    static constexpr bool isComplexValued = GetGetReturnType<T>::isComplex;
    using floatType = typename RadialProjectionResult<sType>::floatType;
    using resultType = RadialProjectionResult<sType>;

    static constexpr size_t NDim = GetNDim::get<T>();

    RadialProjector(const T &instance, SpaceStateType spaceType, device::memory::host_ptr<MemoryToolBox<NDim>> pToolBox,
                    bool pUseCentralBinValues)
        : mSpaceType(spaceType), mInstance(instance), mToolBox(pToolBox), mUseBinCentralValues(pUseCentralBinValues),
          mLayout((mSpaceType == SpaceStateType::Fourier) ? mToolBox->mLayouts.getFourierSpaceLayout()
                                                          : mToolBox->mLayouts.getConfigSpaceLayout())
    {
    }

    /** @brief Measure the projection of your n-dimensional lattice onto the
     *  radial direction. Only supports linear binning, so if you want to
     *  transform to e.g. logarithmic binning, you will simply have a pass
     *  a large number of linear bins here, and do the logarithmic binning yourself.
     *  Default is nGrid ^ (nDimensions - 1).
     */
    RadialProjectionResult<sType> measure(ptrdiff_t nLinearBins = -1, sType customRange = -1, bool excludeOrigin = true)
    {
      if (nLinearBins < 0) {
        ptrdiff_t nGrid = mLayout.getGlobalSizes()[0];
        nLinearBins = std::pow(nGrid, std::max((size_t)1, NDim - 1));
      }

      RadialProjectionResult<sType> baseWorkSpace(nLinearBins, mUseBinCentralValues,
                                                  mSpaceType == SpaceStateType::Fourier);

      sType minValue = excludeOrigin ? 1.0 : 0.0;

      RadialProjectionResult<sType> myResult =
          mSpaceType == SpaceStateType::Configuration
              ? computeConfigurationSpace(makeBinComputer(nLinearBins, minValue, customRange), baseWorkSpace,
                                          excludeOrigin)
              : computeFourierSpace(makeBinComputer(nLinearBins, minValue, customRange), baseWorkSpace, excludeOrigin);

      myResult.finalize(mToolBox->mGroup.getBaseComm());

      return myResult;
    }

  private:
    const SpaceStateType mSpaceType;
    T mInstance;
    device::memory::host_ptr<MemoryToolBox<NDim>> mToolBox;
    bool mUseBinCentralValues;
    LayoutStruct<NDim> mLayout;

  public:
    template <typename BINCOMPUTETYPE>
    RadialProjectionResult<sType> computeConfigurationSpace(BINCOMPUTETYPE binComputer,
                                                            RadialProjectionResult<sType> baseWorkSpace,
                                                            bool excludeOrigin)
    {
      confirmGetterSpace();

#ifdef DEVICE_KOKKOS
      // --- TeamPolicy path: team-local scratch bins ---

      const auto memSizes = mLayout.getSizesInMemory();
      const auto nGhosts = static_cast<device::Idx>(mLayout.getNGhosts());
      ptrdiff_t totalSites = 1;
      for (size_t d = 0; d < NDim; ++d)
        totalSites *= memSizes[d];

      const ptrdiff_t nBins = baseWorkSpace.getNBins();

      const size_t scratchBytes = ScratchBinLayout<sType>::bytesNeeded(nBins);
      constexpr size_t kMaxLevel0Bytes = 48u * 1024u;
      const int scratchLevel = (scratchBytes <= kMaxLevel0Bytes) ? 0 : 1;

      using ExecutionSpace = device_kokkos::DefaultExecutionSpace;
      using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
      using MemberType = typename TeamPolicy::member_type;

      constexpr ptrdiff_t kTargetSitesPerTeam = 1024;
      ptrdiff_t leagueSize = (totalSites + kTargetSitesPerTeam - 1) / kTargetSitesPerTeam;
      if (leagueSize < 1) leagueSize = 1;

      TeamPolicy teamPolicy(leagueSize, Kokkos::AUTO);
      teamPolicy.set_scratch_size(scratchLevel, Kokkos::PerTeam(scratchBytes));

      const auto gVAvg = baseWorkSpace.valuesQuantity().averagesDevice();
      const auto gVVar = baseWorkSpace.valuesQuantity().variancesDevice();
      const auto gVMin = baseWorkSpace.valuesQuantity().minsDevice();
      const auto gVMax = baseWorkSpace.valuesQuantity().maxsDevice();
      const auto gBAvg = baseWorkSpace.binBoundsQuantity().averagesDevice();
      const auto gBVar = baseWorkSpace.binBoundsQuantity().variancesDevice();
      const auto gBMin = baseWorkSpace.binBoundsQuantity().minsDevice();
      const auto gBMax = baseWorkSpace.binBoundsQuantity().maxsDevice();
      const auto gMult = baseWorkSpace.multiplicitiesDevice();

      Kokkos::parallel_for("RadialProjectorConfiguration", teamPolicy,
          DEVICE_CLASS_LAMBDA(const MemberType& team) {

          ScratchBinLayout<sType> scratch(team.team_scratch(scratchLevel), nBins);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nBins), [&](ptrdiff_t b) {
            scratch.init(b);
          });
          team.team_barrier();

          const ptrdiff_t sitesPerTeam = (totalSites + leagueSize - 1) / leagueSize;
          const ptrdiff_t teamStart = team.league_rank() * sitesPerTeam;
          const ptrdiff_t teamEnd = (teamStart + sitesPerTeam < totalSites)
                                        ? teamStart + sitesPerTeam
                                        : totalSites;

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, teamEnd - teamStart),
              [&](ptrdiff_t localIdx) {
            const ptrdiff_t flatIdx = teamStart + localIdx;
            const auto idx = flatToMemoryIndex<NDim>(flatIdx, memSizes, nGhosts);

            device::IdxArray<NDim> global_coords;
            device::apply([&](auto&&... args) {
              mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...);
            }, idx);

            // Origin check (no Hermitian check in configuration space)
            bool isAtOrigin = true;
            for (auto&& it : global_coords)
              if (it != 0) { isAtOrigin = false; break; }
            if (excludeOrigin && isAtOrigin) return;

            sType r{};
            for (size_t i = 0; i < NDim; ++i)
              r += global_coords[i] * global_coords[i];
            r = device::sqrt(r);

            const ptrdiff_t bin = binComputer(r);

            const sType value = device::apply([&](auto&&... args) {
              return DoEval::eval(mInstance, args...);
            }, idx);

            scratch.accumulate(bin, value, r, sType(1));  // weight = 1 always
          });
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nBins), [&](ptrdiff_t b) {
            scratch.mergeTo(gVAvg, gVVar, gVMin, gVMax,
                            gBAvg, gBVar, gBMin, gBMax,
                            gMult, b);
          });
      });

#else
      // --- Original path for DEVICE_STD ---
      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::IdxArray<NDim> global_coords;
        device::apply([&](auto &&...args) { mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...); },
                      idx);

        bool isAtOrigin = true;
        for (auto &&it : global_coords)
          if (it != 0) {
            isAtOrigin = false;
            break;
          }
        if (excludeOrigin && isAtOrigin) [[unlikely]]
          return;

        sType r{};
        for (size_t i = 0; i < NDim; ++i)
          r += global_coords[i] * global_coords[i];
        r = device::sqrt(r);

        const ptrdiff_t bin = binComputer(r);

        device::apply([&](auto &&...args) { baseWorkSpace.add_device(bin, DoEval::eval(mInstance, args...), r, 1.); },
                      idx);
      };

      device::iteration::foreach ("RadialProjectorConfiguration", mLayout, functor);
#endif

      baseWorkSpace.pull();
      binComputer.setCentralBinBounds(baseWorkSpace.getCentralBinBounds());
      return baseWorkSpace;
    }

  public:
    template <typename BINCOMPUTETYPE>
    RadialProjectionResult<sType> computeFourierSpace(BINCOMPUTETYPE binComputer,
                                                      RadialProjectionResult<sType> baseWorkSpace, bool excludeOrigin)
    {
      confirmGetterSpace();

#ifdef DEVICE_KOKKOS
      // --- TeamPolicy path: team-local scratch bins ---

      // 1. Compute flat iteration range
      const auto memSizes = mLayout.getSizesInMemory();
      const auto nGhosts = static_cast<device::Idx>(mLayout.getNGhosts());
      ptrdiff_t totalSites = 1;
      for (size_t d = 0; d < NDim; ++d)
        totalSites *= memSizes[d];

      const ptrdiff_t nBins = baseWorkSpace.getNBins();

      // 2. Determine scratch level: Level 0 (shared memory) if bins fit, else Level 1 (global memory)
      const size_t scratchBytes = ScratchBinLayout<sType>::bytesNeeded(nBins);
      constexpr size_t kMaxLevel0Bytes = 48u * 1024u;
      const int scratchLevel = (scratchBytes <= kMaxLevel0Bytes) ? 0 : 1;

      // 3. Configure TeamPolicy
      using ExecutionSpace = device_kokkos::DefaultExecutionSpace;
      using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
      using MemberType = typename TeamPolicy::member_type;

      constexpr ptrdiff_t kTargetSitesPerTeam = 1024;
      ptrdiff_t leagueSize = (totalSites + kTargetSitesPerTeam - 1) / kTargetSitesPerTeam;
      if (leagueSize < 1) leagueSize = 1;

      TeamPolicy teamPolicy(leagueSize, Kokkos::AUTO);
      teamPolicy.set_scratch_size(scratchLevel, Kokkos::PerTeam(scratchBytes));

      // 4. Extract device view handles for the merge step
      const auto gVAvg = baseWorkSpace.valuesQuantity().averagesDevice();
      const auto gVVar = baseWorkSpace.valuesQuantity().variancesDevice();
      const auto gVMin = baseWorkSpace.valuesQuantity().minsDevice();
      const auto gVMax = baseWorkSpace.valuesQuantity().maxsDevice();
      const auto gBAvg = baseWorkSpace.binBoundsQuantity().averagesDevice();
      const auto gBVar = baseWorkSpace.binBoundsQuantity().variancesDevice();
      const auto gBMin = baseWorkSpace.binBoundsQuantity().minsDevice();
      const auto gBMax = baseWorkSpace.binBoundsQuantity().maxsDevice();
      const auto gMult = baseWorkSpace.multiplicitiesDevice();

      Kokkos::parallel_for("RadialProjectorFourier", teamPolicy,
          DEVICE_CLASS_LAMBDA(const MemberType& team) {

          // --- Step A: Allocate scratch bins ---
          ScratchBinLayout<sType> scratch(team.team_scratch(scratchLevel), nBins);

          // --- Step B: Initialize scratch to identity values ---
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nBins), [&](ptrdiff_t b) {
            scratch.init(b);
          });
          team.team_barrier();

          // --- Step C: Accumulate lattice sites into scratch bins ---
          const ptrdiff_t sitesPerTeam = (totalSites + leagueSize - 1) / leagueSize;
          const ptrdiff_t teamStart = team.league_rank() * sitesPerTeam;
          const ptrdiff_t teamEnd = (teamStart + sitesPerTeam < totalSites)
                                        ? teamStart + sitesPerTeam
                                        : totalSites;

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, teamEnd - teamStart),
              [&](ptrdiff_t localIdx) {
            const ptrdiff_t flatIdx = teamStart + localIdx;

            // Flat → ND memory indices (with ghost offset)
            const auto idx = flatToMemoryIndex<NDim>(flatIdx, memSizes, nGhosts);

            // Get global spatial coordinates
            device::IdxArray<NDim> global_coords;
            device::apply([&](auto&&... args) {
              mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...);
            }, idx);

            // Origin check
            bool isAtOrigin = true;
            for (auto&& it : global_coords)
              if (it != 0) { isAtOrigin = false; break; }
            if (excludeOrigin && isAtOrigin) return;

            // Hermitian partner check
            const HermitianRedundancy quality = mLayout.getHermitianPartners().qualify(global_coords);
            if (quality == HermitianRedundancy::negativePartner) return;

            // Compute radius
            sType r{};
            for (size_t i = 0; i < NDim; ++i)
              r += global_coords[i] * global_coords[i];
            r = device::sqrt(r);

            // Map radius to bin
            const ptrdiff_t bin = binComputer(r);

            // Weight: 0.5 for real-valued entries, 1.0 otherwise
            const sType weight = (quality == HermitianRedundancy::realValued) ? sType(0.5) : sType(1);

            // Evaluate the field expression at this site
            const sType value = device::apply([&](auto&&... args) {
              return DoEval::eval(mInstance, args...);
            }, idx);

            // Accumulate into team-local scratch
            scratch.accumulate(bin, value, r, weight);
          });
          team.team_barrier();

          // --- Step D: Merge scratch into global arrays ---
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nBins), [&](ptrdiff_t b) {
            scratch.mergeTo(gVAvg, gVVar, gVMin, gVMax,
                            gBAvg, gBVar, gBMin, gBMax,
                            gMult, b);
          });
      });

#else
      // --- Original path for DEVICE_STD (serial, no Kokkos) ---
      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        device::IdxArray<NDim> global_coords;
        device::apply([&](auto &&...args) { mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...); },
                      idx);

        bool isAtOrigin = true;
        for (auto &&it : global_coords)
          if (it != 0) {
            isAtOrigin = false;
            break;
          }
        if (excludeOrigin && isAtOrigin) [[unlikely]]
          return;

        const HermitianRedundancy quality = mLayout.getHermitianPartners().qualify(global_coords);
        if (quality != HermitianRedundancy::negativePartner) {
          sType r{};
          for (size_t i = 0; i < NDim; ++i)
            r += global_coords[i] * global_coords[i];
          r = device::sqrt(r);

          const ptrdiff_t bin = binComputer(r);
          floatType weight = quality == HermitianRedundancy::realValued ? 0.5 : 1;

          device::apply(
              [&](auto &&...args) { baseWorkSpace.add_device(bin, DoEval::eval(mInstance, args...), r, weight); }, idx);
        }
      };
      device::iteration::foreach ("RadialProjectorFourier", mLayout, functor);
#endif

      baseWorkSpace.pull();
      binComputer.setCentralBinBounds(baseWorkSpace.getCentralBinBounds());

      return baseWorkSpace;
    }

  private:
    void confirmGetterSpace()
    {
      ConfirmSpace::apply(mInstance, mLayout, mSpaceType);
      GhostsHunter::apply(mInstance);
    }

    /** @brief Creates the lambda that maps the IterationCoordinates to a bin. */
    inline auto makeBinComputer(ptrdiff_t nLinearBins, sType minValue, sType customRange = -1)
    {
      auto rMax = customRange < 0 ? mLayout.getMaxRadius() : customRange;

      return RadialBinComputer(minValue, rMax, nLinearBins);
    }
  };

  template <size_t NDim, typename T>
  RadialProjector<T> projectRadially(T instance, SpaceStateType spaceType,
                                     device::memory::host_ptr<MemoryToolBox<NDim>> pToolBox, bool useBinCentralValues)
  {
    return RadialProjector<T>(instance, spaceType, pToolBox, useBinCentralValues);
  }

  template <typename T> RadialProjector<T> projectRadially(T instance, bool useBinCentralValues = false)
  {
    return projectRadially(instance, SpaceStateType::Configuration, GetToolBox::get(instance), useBinCentralValues);
  }

  template <typename T> RadialProjector<T> projectRadiallyFourier(T instance, bool useBinCentralValues = false)
  {
    return projectRadially(instance, SpaceStateType::Fourier, GetToolBox::get(instance), useBinCentralValues);
  }

} // namespace TempLat

#endif
