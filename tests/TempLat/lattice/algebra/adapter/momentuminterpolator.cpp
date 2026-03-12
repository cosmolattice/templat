/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/adapter/momentuminterpolator.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"

#include <cmath>
#include <vector>

namespace TempLat
{
  template <size_t NDim> struct MomentumInterpolatorTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> inline void MomentumInterpolatorTester<NDim>::Test(TDDAssertion &tdd)
  {
    const ptrdiff_t nGrid = 16, nGhost = 0;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // Build a known function: ps(k) = k^2 on a fine grid covering the range of
    // possible momenta. The maximum dimensionless wavenumber magnitude is
    // sqrt(NDim) * (nGrid/2), so with kIR=1 we need to cover [0, that].
    const double kIR = 1.0;
    const double kMaxSpline = std::sqrt(static_cast<double>(NDim)) * (nGrid / 2) + 1.0;
    const int nPoints = 200;
    std::vector<double> kVec(nPoints), psVec(nPoints);
    for (int i = 0; i < nPoints; ++i) {
      kVec[i] = kMaxSpline * i / (nPoints - 1);
      psVec[i] = kVec[i] * kVec[i]; // ps = k^2
    }

    MomentumInterpolator<double, NDim> interp(kVec, psVec, toolBox, kIR);

    // Assign to a field in Fourier space
    Field<double, NDim> phi("phi", toolBox);
    phi.inFourierSpace() = interp;

    // Read back and verify
    auto phi_view = phi.inFourierSpace().getLocalNDHostView();
    auto layout = toolBox->mLayouts.getFourierSpaceLayout();

    bool correct = true;

    NDLoop<NDim>(phi_view, [&](const auto &...indices) {
      // Compute the expected value: |ntilde| * kIR, then ps(k) = k^2
      device::IdxArray<NDim> global_idx;
      layout.putSpatialLocationFromMemoryIndexInto(global_idx, indices...);

      double norm2 = 0.0;
      for (size_t d = 0; d < NDim; ++d) {
        double n = static_cast<double>(global_idx[d]);
        norm2 += n * n;
      }
      double k = std::sqrt(norm2) * kIR;

      // The interpolator clamps to [kMin, kMax] and evaluates the spline.
      // For ps = k^2 with a fine cubic spline, the result should be very close to k^2.
      double expected = k * k;
      double actual = phi_view(indices...).real();

      // Cubic spline interpolation of k^2 should be nearly exact
      bool this_correct = std::fabs(actual - expected) < 0.1 * (1.0 + expected);
      correct &= this_correct;
    });

    tdd.verify(correct, "MomentumInterpolator correctly maps Fourier sites through spline");

    // Test with a different kIR
    const double kIR2 = 0.5;
    // Need a spline covering [0, kMaxSpline * kIR2]
    const double kMaxSpline2 = std::sqrt(static_cast<double>(NDim)) * (nGrid / 2) * kIR2 + 1.0;
    std::vector<double> kVec2(nPoints), psVec2(nPoints);
    for (int i = 0; i < nPoints; ++i) {
      kVec2[i] = kMaxSpline2 * i / (nPoints - 1);
      psVec2[i] = std::sin(kVec2[i]); // ps = sin(k)
    }

    MomentumInterpolator<double, NDim> interp2(kVec2, psVec2, toolBox, kIR2);

    Field<double, NDim> phi2("phi2", toolBox);
    phi2.inFourierSpace() = interp2;

    auto phi2_view = phi2.inFourierSpace().getLocalNDHostView();
    bool correct2 = true;

    NDLoop<NDim>(phi2_view, [&](const auto &...indices) {
      device::IdxArray<NDim> global_idx;
      layout.putSpatialLocationFromMemoryIndexInto(global_idx, indices...);

      double norm2 = 0.0;
      for (size_t d = 0; d < NDim; ++d) {
        double n = static_cast<double>(global_idx[d]);
        norm2 += n * n;
      }
      double k = std::sqrt(norm2) * kIR2;

      // Clamp to spline range
      double kClamped = std::max(kVec2.front(), std::min(kVec2.back(), k));
      double expected = std::sin(kClamped);
      double actual = phi2_view(indices...).real();

      bool this_correct = std::fabs(actual - expected) < 0.01 * (1.0 + std::fabs(expected));
      correct2 &= this_correct;
    });

    tdd.verify(correct2, "MomentumInterpolator with kIR=0.5 and sin(k) spline");

    // Test clamping: use a spline with narrow range [1, 3], so sites with k outside get clamped
    std::vector<double> kNarrow = {1.0, 1.5, 2.0, 2.5, 3.0};
    std::vector<double> psNarrow = {10.0, 12.25, 14.0, 16.25, 19.0};

    MomentumInterpolator<double, NDim> interpNarrow(kNarrow, psNarrow, toolBox, kIR);

    Field<double, NDim> phi3("phi3", toolBox);
    phi3.inFourierSpace() = interpNarrow;

    auto phi3_view = phi3.inFourierSpace().getLocalNDHostView();
    bool correct3 = true;

    // Check the k=0 site (all indices zero in Fourier space). It should clamp to kMin=1.0
    // and return ps(1.0) = 10.0.
    // The zero-frequency site is at the origin of the local Fourier layout.
    // For a simple check, verify that all values are within [ps(kMin), ps(kMax)] approximately.
    NDLoop<NDim>(phi3_view, [&](const auto &...indices) {
      double actual = phi3_view(indices...).real();
      // All values should be >= ps(kMin) - tolerance and <= ps(kMax) + tolerance
      // since clamp brings out-of-range k to the boundary
      bool this_correct = (actual >= 9.5) && (actual <= 19.5);
      correct3 &= this_correct;
    });

    tdd.verify(correct3, "MomentumInterpolator clamping keeps values in spline range");
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MomentumInterpolatorTester<2>> test2;
  TempLat::TDDContainer<TempLat::MomentumInterpolatorTester<3>> test3;
} // namespace
