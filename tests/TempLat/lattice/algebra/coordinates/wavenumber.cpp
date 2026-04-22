/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler,  Year: 2025

#include "TempLat/lattice/algebra/coordinates/wavenumber.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/util/log/saycomplete.h"
#include "TempLat/util/ndloop.h"

namespace TempLat
{

  template <typename T, size_t NDim> struct WaveNumberTester {
    static void Test(TDDAssertion &tdd);
  };

  template <typename T, size_t NDim> inline void WaveNumberTester<T, NDim>::Test(TDDAssertion &tdd)
  {
    const device::Idx nGrid = 16, nGhost = 0;

    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // Create fields for each wavenumber component
    std::vector<Field<T, NDim>> phi_components;
    phi_components.reserve(NDim);
    for (size_t d = 0; d < NDim; ++d) {
      phi_components.emplace_back("phi_" + std::to_string(d), toolBox);
    }

    Field<T, NDim> phinorm("phinorm", toolBox);
    Field<T, NDim> phinorm2("phinorm2", toolBox);

    WaveNumber<NDim> k(toolBox);

    // Assign wavenumber components to fields
    constexpr_for<0, NDim>([&](auto d) { phi_components[d].inFourierSpace() = k(d + Tag<1>()); });
    phinorm.inFourierSpace() = k.norm();
    phinorm2.inFourierSpace() = k.norm2();

    // Get host views for all component fields
    std::vector<decltype(phi_components[0].inFourierSpace().getLocalNDHostView())> phi_views;
    phi_views.reserve(NDim);
    for (size_t d = 0; d < NDim; ++d) {
      phi_views.push_back(phi_components[d].inFourierSpace().getLocalNDHostView());
    }
    auto phinorm_view = phinorm.inFourierSpace().getLocalNDHostView();
    auto phinorm2_view = phinorm2.inFourierSpace().getLocalNDHostView();

    // Get Fourier space layout for MPI support (handles transposition)
    auto layout = toolBox->mLayouts.getFourierSpaceLayout();

    // Check that the fourier coordinate is correct
    bool correct = true;

    NDLoop<NDim>(phinorm_view, [&](const auto &...indices) {
      device::IdxArray<NDim> local_idx{indices...};

      // Use layout to compute global spatial coordinates from local memory indices
      // This correctly handles transposition in MPI distributed FFTs
      device::IdxArray<NDim> global_idx;
      layout.putSpatialLocationFromMemoryIndexInto(global_idx, indices...);

      bool this_correct = true;

      // Check each wavenumber component - the value should match the global spatial coordinate
      T norm2 = 0.0;
      for (size_t d = 0; d < NDim; ++d) {
        const T expected_val = static_cast<T>(global_idx[d]);
        norm2 += expected_val * expected_val;

        this_correct &= AlmostEqual(T(phi_views[d](indices...).real()), T(expected_val));
        this_correct &= AlmostEqual(T(phi_views[d](indices...).imag()), T(0.));
      }

      std::stringstream ss;
      // Check norm2 and norm
      this_correct &= AlmostEqual(T(phinorm2_view(indices...).real()), T(norm2));
      ss << "Fail after 1, should: " << norm2 << " is " << phinorm2_view(indices...).real() << "\n";
      this_correct &= AlmostEqual(T(phinorm2_view(indices...).imag()), T(0.));
      ss << "Fail after 2, should: " << 0. << " is " << phinorm2_view(indices...).imag() << "\n";
      this_correct &= AlmostEqual(T(phinorm_view(indices...).real()), T(sqrt(norm2)));
      ss << "Fail after 3, should: " << sqrt(norm2) << " is " << phinorm_view(indices...).real() << "\n";
      this_correct &= AlmostEqual(T(phinorm_view(indices...).imag()), T(0.));
      ss << "Fail after 4, should: " << 0. << " is " << phinorm_view(indices...).imag() << "\n";

      correct &= this_correct;

      if (!this_correct) {
        ss << "Error at local (";
        for (size_t d = 0; d < NDim; ++d) {
          ss << local_idx[d];
          if (d < NDim - 1) ss << ", ";
        }
        ss << "), global (";
        for (size_t d = 0; d < NDim; ++d) {
          ss << global_idx[d];
          if (d < NDim - 1) ss << ", ";
        }
        ss << "): ";
        for (size_t d = 0; d < NDim; ++d) {
          ss << "phi[" << d << "] = " << phi_views[d](indices...);
          if (d < NDim - 1) ss << ", ";
        }
        ss << ", phinorm = " << phinorm_view(indices...) << ", phinorm2 = " << phinorm2_view(indices...) << "\n";
        sayMPI << ss.str();
      }
    });

    tdd.verify(correct);
  }

} // namespace TempLat

namespace
{
#ifndef HAVE_FFTFLOAT
  TempLat::TDDContainer<TempLat::WaveNumberTester<double, 2>> test2;
  TempLat::TDDContainer<TempLat::WaveNumberTester<double, 3>> test3;
  TempLat::TDDContainer<TempLat::WaveNumberTester<double, 4>> test4;
  TempLat::TDDContainer<TempLat::WaveNumberTester<double, 5>> test5;
#else
  TempLat::TDDContainer<TempLat::WaveNumberTester<float, 2>> test6;
  TempLat::TDDContainer<TempLat::WaveNumberTester<float, 3>> test7;
  TempLat::TDDContainer<TempLat::WaveNumberTester<float, 4>> test8;
  TempLat::TDDContainer<TempLat::WaveNumberTester<float, 5>> test9;
#endif

} // namespace
