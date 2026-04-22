/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/coordinates/momentummultiplicity.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/util/log/saycomplete.h"
#include "TempLat/util/ndloop.h"

namespace TempLat
{

  template <size_t NDim> struct MomentumMultiplicityTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> inline void MomentumMultiplicityTester<NDim>::Test(TDDAssertion &tdd)
  {
    const device::Idx nGrid = 16, nGhost = 0;
    const device::Idx Nh = nGrid / 2;

    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // --- Test 1: Bin counts sum to total non-zero Fourier modes ---
    auto binP = getTypeIBinCounts<NDim>(nGrid);

    // Total non-zero modes in the full grid = N^D - 1
    device::Idx totalModes = 1;
    for (size_t d = 0; d < NDim; ++d)
      totalModes *= nGrid;
    totalModes -= 1; // exclude the zero mode

    device::Idx binSum = 0;
    for (auto &b : binP)
      binSum += b;

    tdd.verify(binSum == totalModes);

    // --- Test 2: Field assignment via expression templates ---
    MomentumMultiplicity<double, NDim> mult(toolBox);
    auto layout = toolBox->mLayouts.getFourierSpaceLayout();

    Field<double, NDim> phi("phi_mult", toolBox);
    phi.inFourierSpace() = mult;

    auto phi_view = phi.inFourierSpace().getLocalNDHostView();

    bool fieldCorrect = true;

    NDLoop<NDim>(phi_view, [&](const auto &...indices) {
      const double val = phi_view(indices...).real();

      device::IdxArray<NDim> global_idx;
      layout.putSpatialLocationFromMemoryIndexInto(global_idx, indices...);

      size_t n2 = 0;
      bool allZero = true;
      bool allHalf = true;
      for (size_t d = 0; d < NDim; ++d) {
        const size_t folded = static_cast<size_t>(std::abs(global_idx[d]));
        n2 += folded * folded;
        if (folded != 0) allZero = false;
        if (folded != static_cast<size_t>(Nh)) allHalf = false;
      }

      double expected;
      if (allZero) {
        expected = 1.0;
      } else if (allHalf) {
        expected = 1.0 / double(1 << NDim);
      } else {
        const int bin = static_cast<int>(std::sqrt(double(n2)) + 0.5) - 1;
        expected = 1.0 / double(binP[bin]);
      }

      bool ok = AlmostEqual(val, expected) && AlmostEqual(phi_view(indices...).imag(), 0.0);
      if (!ok) {
        std::stringstream ss;
        ss << "Mismatch: field=" << val << " expected=" << expected << " imag=" << phi_view(indices...).imag()
           << " n2=" << n2;
        ss << " coords=(";
        for (size_t d = 0; d < NDim; ++d) {
          ss << global_idx[d];
          if (d < NDim - 1) ss << ",";
        }
        ss << ")";
        sayMPI << ss.str();
      }
      fieldCorrect &= ok;
    });

    tdd.verify(fieldCorrect);
  }

} // namespace TempLat

namespace
{
  // NDim=1 skipped due to FFT/MPI issue (same as wavenumber test)
  TempLat::TDDContainer<TempLat::MomentumMultiplicityTester<2>> test2;
  TempLat::TDDContainer<TempLat::MomentumMultiplicityTester<3>> test3;
  TempLat::TDDContainer<TempLat::MomentumMultiplicityTester<4>> test4;
} // namespace
