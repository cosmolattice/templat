#ifndef TEMPLAT_LATTICE_ALGEBRA_COORDINATES_MULTIPLICITY_H
#define TEMPLAT_LATTICE_ALGEBRA_COORDINATES_MULTIPLICITY_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Nicolas Loayza, Franz R. Sattler,  Year: 2025

#include <vector>
#include <cmath>
#include <cstddef>
#include <utility>

#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/spacestateinterface.h"
#include "TempLat/util/exception.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/util/constexpr_for.h"

namespace TempLat
{
  MakeException(MomentumMultiplicityWrongSpaceConfirmation);

  namespace detail
  {
    /**
     * @brief Recursively iterates over the positive octant [0, Nh]^NDim, accumulating
     *        bin counts with symmetry factors.
     */
    template <size_t Dim, size_t NDim, typename T>
    void iterateOctant(const ptrdiff_t Nh, const ptrdiff_t N, std::vector<int> &binP, ptrdiff_t n2, int nEdge)
    {
      for (ptrdiff_t a = 0; a <= Nh; ++a) {
        const ptrdiff_t newN2 = n2 + a * a;
        const int newNEdge = nEdge + ((a == 0 || a == Nh) ? 1 : 0);

        if constexpr (Dim + 1 < NDim) {
          iterateOctant<Dim + 1, NDim, T>(Nh, N, binP, newN2, newNEdge);
        } else {
          // Innermost dimension: compute bin and accumulate
          if (newN2 == 0) continue; // skip the all-zero mode

          const int bin = static_cast<int>(std::sqrt(T(newN2)) + T(0.5)) - 1;
          binP[bin] += ((1 << NDim) >> newNEdge);
        }
      }
    }
  } // namespace detail

  template <size_t NDim> std::vector<int> getTypeIBinCounts(const ptrdiff_t N)
  {
    const auto Nh = N / 2;

    // Max bin index: floor(sqrt(NDim) * Nh) + 1
    const int maxBin = static_cast<int>(std::sqrt(double(NDim)) * double(Nh)) + 1;
    std::vector<int> binP(maxBin, 0);

    detail::iterateOctant<0, NDim, double>(Nh, N, binP, 0, 0);

    return binP;
  }

  /**
   * @brief Computes, for each Fourier-space site, the multiplicity of modes sharing
   *        the same rounded momentum magnitude bin. Returns 1/multiplicity on eval.
   *        Generalized to arbitrary NDim.
   *
   * Unit test: ctest -R test-momentummultiplicity
   */
  template <typename T, size_t NDim> class MomentumMultiplicity
  {
  public:
    MomentumMultiplicity(device::memory::host_ptr<MemoryToolBox<NDim>> toolBox)
        : mnGrid(toolBox->mNGridPointsVec[0]), mLayout(toolBox->mLayouts.getFourierSpaceLayout())
    {
      const auto binP = getTypeIBinCounts<NDim>(mnGrid);
      // Allocate device-accessible view and copy from host
      mmulti = device::memory::NDView<int, 1>("MomentumMultiplicity_bins", binP.size());
      device::memory::copyHostToDevice(binP.data(), mmulti);
    }

    DEVICE_FUNCTION
    MomentumMultiplicity(const MomentumMultiplicity &) = default;

    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_INLINE_FUNCTION T eval(const IDX &...idx) const
    {
      // Convert memory indices to spatial coordinates
      device::IdxArray<NDim> spatial;
      mLayout.putSpatialLocationFromMemoryIndexInto(spatial, idx...);

      const ptrdiff_t half = mnGrid / 2;

      // Fold spatial coordinates (which may be negative) to [0, Nh]
      size_t n2 = 0;
      bool allZero = true;
      bool allHalf = true;
      for (size_t d = 0; d < NDim; ++d) {
        const size_t folded = static_cast<size_t>((spatial[d] >= 0) ? spatial[d] : -spatial[d]);
        n2 += folded * folded;
        if (folded != 0) allZero = false;
        if (folded != static_cast<size_t>(half)) allHalf = false;
      }

      T pair = T(1);

      if (!allZero) {
        const int bin = static_cast<int>(std::sqrt(T(n2)) + T(0.5)) - 1;
        pair = T(mmulti(bin));
      }

      if (allHalf) {
        pair = T(1 << NDim);
      }

      return T(1) / pair;
    }

    static std::string toString() { return "MomentumMultiplicity"; }

    void confirmSpace(ptrdiff_t i, const LayoutStruct<NDim> &newLayout, const SpaceStateType &spaceType) const
    {
      switch (spaceType) {
      case SpaceStateType::Configuration:
        throw MomentumMultiplicityWrongSpaceConfirmation(
            "MomentumMultiplicity explicitly only can be used in fourier space. Abort.");
        break;
      case SpaceStateType::Fourier:
      default:
        break;
      }
    }

  private:
    ptrdiff_t mnGrid;
    LayoutStruct<NDim> mLayout;
    device::memory::NDView<int, 1> mmulti;
  };

} // namespace TempLat

#endif
