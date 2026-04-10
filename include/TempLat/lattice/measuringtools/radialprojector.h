#ifndef TEMPLAT_LATTICE_MEASUREMENTS_RADIALPROJECTOR_H
#define TEMPLAT_LATTICE_MEASUREMENTS_RADIALPROJECTOR_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019
//            Modified by: Jorge Baeza-Ballesteros, Year: 2024

#include <algorithm>
#include <cstddef>

#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionresult.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/unbinnedradialprojectionresult.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglequantity.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialbincomputer.h"
#include "TempLat/lattice/algebra/operators/squareroot.h"

#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/algebra/spacestateinterface.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/lattice/algebra/helpers/getngrid.h"

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
    RadialProjectionResult<sType> measure(sType maxValue, sType deltakBins = 1, bool excludeOrigin = true)
    {
      sType minValue = !excludeOrigin ? 0.0 : deltakBins >= 1.0 ? 0.5 : 1.0 - deltakBins / 2.;
      ptrdiff_t nLinearBins = ceil( (maxValue - minValue) / deltakBins );

      RadialProjectionResult<sType> baseWorkSpace(nLinearBins, mUseBinCentralValues,
                                                  mSpaceType == SpaceStateType::Fourier);

      RadialProjectionResult<sType> myResult =
          mSpaceType == SpaceStateType::Configuration
              ? computeConfigurationSpace(makeBinComputer(nLinearBins, minValue, maxValue, deltakBins), baseWorkSpace,
                                          excludeOrigin)
              : computeFourierSpace(makeBinComputer(nLinearBins, minValue, maxValue, deltakBins), baseWorkSpace, excludeOrigin);

      myResult.finalize(mToolBox->mGroup.getBaseComm());

      return myResult;
    }

    UnbinnedRadialProjectionResult<sType> measureUnbinned(ptrdiff_t N, bool excludeOrigin = true) {

      ptrdiff_t nLinearBins = NDim * pow<2>(N) / 4 + 1;

      UnbinnedRadialProjectionResult<sType> baseWorkSpace(nLinearBins,  mSpaceType == SpaceStateType::Fourier);


      UnbinnedRadialProjectionResult<sType> myResult = mSpaceType == SpaceStateType::Configuration ?
      computeConfigurationSpaceUnbinned(baseWorkSpace, excludeOrigin) :
      computeFourierSpaceUnbinned(baseWorkSpace, excludeOrigin);


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

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        // Get the global coordinates of this index.
        device::IdxArray<NDim> global_coords;
        device::apply([&](auto &&...args) { mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...); },
                      idx);

        // Check if we are at the origin.
        bool isAtOrigin = true;
        for (auto &&it : global_coords) {
          if (it != 0) {
            isAtOrigin = false;
            break;
          }
        }
        if (excludeOrigin && isAtOrigin) [[unlikely]]
          return;

        // get the radius
        sType r{};
        for (size_t i = 0; i < NDim; ++i)
          r += global_coords[i] * global_coords[i];
        r = device::sqrt(r);

        // Map the radius to a bin
        const ptrdiff_t bin = binComputer(r);

        // Add the bin contribution to the workspace.
        device::apply([&](auto &&...args) { baseWorkSpace.add_device(bin, DoEval::eval(mInstance, args...), r, 1.); },
                      idx);
      };

      device::iteration::foreach ("RadialProjectorConfiguration", mLayout, functor);

      baseWorkSpace.pull();
      binComputer.setCentralBinBounds(baseWorkSpace.getCentralBinBounds());
      return baseWorkSpace;
    }

    template <typename BINCOMPUTETYPE>
    RadialProjectionResult<sType> computeFourierSpace(BINCOMPUTETYPE binComputer,
                                                      RadialProjectionResult<sType> baseWorkSpace, bool excludeOrigin)
    {
      confirmGetterSpace();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        // Get the global coordinates of this index.
        device::IdxArray<NDim> global_coords;
        device::apply([&](auto &&...args) { mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...); },
                      idx);

        // Check if we are at the origin.
        bool isAtOrigin = true;
        for (auto &&it : global_coords) {
          if (it != 0) {
            isAtOrigin = false;
            break;
          }
        }
        if (excludeOrigin && isAtOrigin) [[unlikely]]
          return;

        const HermitianRedundancy quality = mLayout.getHermitianPartners().qualify(global_coords);
        if (quality != HermitianRedundancy::negativePartner) {
          // get the radius
          sType r{};
          for (size_t i = 0; i < NDim; ++i)
            r += global_coords[i] * global_coords[i];
          r = device::sqrt(r);

          // Map the radius to a bin
          const ptrdiff_t bin = binComputer(r);

          // don't over-weight the real-valued entries: only one float value, only half the weight.
          floatType weight = quality == HermitianRedundancy::realValued ? 0.5 : 1;

          // Add the bin contribution to the workspace.
          device::apply(
              [&](auto &&...args) { baseWorkSpace.add_device(bin, DoEval::eval(mInstance, args...), r, weight); }, idx);
        }
      };
      device::iteration::foreach ("RadialProjectorFourier", mLayout, functor);
      baseWorkSpace.pull();
      binComputer.setCentralBinBounds(baseWorkSpace.getCentralBinBounds());

      return baseWorkSpace;
    }


    UnbinnedRadialProjectionResult<sType> computeConfigurationSpaceUnbinned(UnbinnedRadialProjectionResult<sType> baseWorkSpace, bool excludeOrigin)
    {
      confirmGetterSpace();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        // Get the global coordinates of this index.
        device::IdxArray<NDim> global_coords;
        device::apply([&](auto &&...args) { mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...); },
                      idx);

        // Check if we are at the origin.
        bool isAtOrigin = true;
        for (auto &&it : global_coords) {
          if (it != 0) {
            isAtOrigin = false;
            break;
          }
        }
        if (excludeOrigin && isAtOrigin) [[unlikely]]
          return;

        // get the radius
        sType r2{};
        for (size_t i = 0; i < NDim; ++i)
          r2 += global_coords[i] * global_coords[i];

        // Add the bin contribution to the workspace.
        device::apply([&](auto &&...args) { baseWorkSpace.add_device(r2, DoEval::eval(mInstance, args...), 1.); },
                      idx);
      };

      device::iteration::foreach ("UnbinnedRadialProjectorConfiguration", mLayout, functor);

      baseWorkSpace.pull();
      return baseWorkSpace;
    }

    UnbinnedRadialProjectionResult<sType> computeFourierSpaceUnbinned(UnbinnedRadialProjectionResult<sType> baseWorkSpace, bool excludeOrigin)
    {
      confirmGetterSpace();

      auto functor = DEVICE_CLASS_LAMBDA(const device::IdxArray<NDim> &idx)
      {
        // Get the global coordinates of this index.
        device::IdxArray<NDim> global_coords;
        device::apply([&](auto &&...args) { mLayout.putSpatialLocationFromMemoryIndexInto(global_coords, args...); },
                      idx);

        // Check if we are at the origin.
        bool isAtOrigin = true;
        for (auto &&it : global_coords) {
          if (it != 0) {
            isAtOrigin = false;
            break;
          }
        }
        if (excludeOrigin && isAtOrigin) [[unlikely]]
          return;

        const HermitianRedundancy quality = mLayout.getHermitianPartners().qualify(global_coords);
        if (quality != HermitianRedundancy::negativePartner) {
          // get the radius
          sType k2{};
          for (size_t i = 0; i < NDim; ++i)
            k2 += global_coords[i] * global_coords[i];

          // don't over-weight the real-valued entries: only one float value, only half the weight.
          floatType weight = quality == HermitianRedundancy::realValued ? 0.5 : 1;

          // Add the bin contribution to the workspace.
          device::apply(
            [&](auto &&...args) { baseWorkSpace.add_device(k2, DoEval::eval(mInstance, args...), weight); }, idx);
        }
      };
      device::iteration::foreach ("RadialProjectorFourier", mLayout, functor);
      baseWorkSpace.pull();

      return baseWorkSpace;
    }

  private:
    void confirmGetterSpace()
    {
      ConfirmSpace::apply(mInstance, mLayout, mSpaceType);
      GhostsHunter::apply(mInstance);
    }

    /** @brief Creates the lambda that maps the IterationCoordinates to a bin. */
    inline auto makeBinComputer(ptrdiff_t nLinearBins, sType minValue, sType maxValue = -1, sType deltakBins = -1)
    {
      auto rMax = maxValue < 0 ? mLayout.getMaxRadius() : maxValue;

      return RadialBinComputer(minValue, rMax, nLinearBins, deltakBins);
    }
  };

  template <size_t NDim, typename T>
  RadialProjector<T> projectRadially(T instance, SpaceStateType spaceType,
                                     device::memory::host_ptr<MemoryToolBox<NDim>> pToolBox, bool useBinCentralValues = false)
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
