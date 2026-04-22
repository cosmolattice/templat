
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#include "TempLat/lattice/memory/memorylayouts/checkerboardlayout.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  template <size_t NDim> struct CheckerboardLayoutTester {
    static void Test(TDDAssertion &tdd);
  };

  /** @brief Helper: compute spatial parity at memory indices, given layout info. */
  template <size_t NDim>
  static device::Idx computeParity(const device::IdxArray<NDim> &memIdx, const device::IdxArray<NDim> &localStarts,
                                   const device::IdxArray<NDim> &transpositionForward, device::Idx nGhosts)
  {
    device::Idx sum = 0;
    for (size_t d = 0; d < NDim; ++d) {
      device::Idx globalDim = transpositionForward[d];
      device::Idx spatialCoord = localStarts[globalDim] + (memIdx[d] - nGhosts);
      sum += spatialCoord;
    }
    return ((sum % 2) + 2) % 2;
  }

  // ---- 3D untransposed test ----
  template <> void CheckerboardLayoutTester<3>::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    const device::Idx nGhosts = 2;
    const device::IdxArray<NDim> gridSizes{{8, 6, 7}}; // 7 is odd in split dim

    LayoutStruct<NDim> layout(gridSizes, nGhosts);
    layout.setLocalStarts({{3, 5, 2}});

    CheckerboardLayout<NDim> even(layout, Parity::Even);
    CheckerboardLayout<NDim> odd(layout, Parity::Odd);

    // Transposition info for parity check
    device::IdxArray<NDim> transForward;
    for (size_t d = 0; d < NDim; ++d)
      transForward[d] = layout.getTranspositionMap_memoryToGlobalSpace().getForward(d);

    // Count sites via reduce
    device::Idx evenCount = 0;
    device::iteration::reduce(
        "CountEven", even, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx, device::Idx &update) { update += 1; },
        evenCount);

    device::Idx oddCount = 0;
    device::iteration::reduce(
        "CountOdd", odd, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx, device::Idx &update) { update += 1; },
        oddCount);

    device::Idx totalSites = 1;
    for (size_t d = 0; d < NDim; ++d)
      totalSites *= gridSizes[d];

    say << "3D untransposed: evenCount=" << evenCount << " oddCount=" << oddCount << " total=" << totalSites << "\n";

    tdd.verify(evenCount + oddCount == totalSites);
    tdd.verify(evenCount > 0);
    tdd.verify(oddCount > 0);

    // Verify parity correctness using a 1D device view
    auto sizesInMem = layout.getSizesInMemory();
    device::Idx flatSize = 1;
    for (size_t d = 0; d < NDim; ++d)
      flatSize *= (sizesInMem[d] + 2 * nGhosts);

    auto parityView = device::memory::NDView<device::Idx, 1>("parityCheck", flatSize);

    // Initialize to -1
    device::iteration::foreach<1>(
        "Init", {0}, {flatSize}, DEVICE_LAMBDA(const device::IdxArray<1> &i) { parityView(i[0]) = -1; });

    // Mark even sites with 0
    const auto localStarts = layout.getLocalStarts();
    device::iteration::foreach (
        "MarkEven", even, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx) {
          // Compute flat index for unique identification
          device::Idx flat = fullIdx[0];
          for (size_t d = 1; d < NDim; ++d)
            flat = flat * (sizesInMem[d] + 2 * nGhosts) + fullIdx[d];
          parityView(flat) = 0;
        });

    // Mark odd sites with 1
    device::iteration::foreach (
        "MarkOdd", odd, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx) {
          device::Idx flat = fullIdx[0];
          for (size_t d = 1; d < NDim; ++d)
            flat = flat * (sizesInMem[d] + 2 * nGhosts) + fullIdx[d];
          parityView(flat) = 1;
        });

    // Copy to host and verify parity + completeness
    std::vector<device::Idx> hostParity(flatSize, -1);
    device::memory::copyDeviceToHost(parityView, hostParity.data());

    device::Idx verifiedEven = 0, verifiedOdd = 0;
    // Iterate over local sites (no ghosts) and check each one
    std::function<void(size_t, device::IdxArray<NDim> &)> iterate;
    iterate = [&](size_t dim, device::IdxArray<NDim> &idx) {
      if (dim == NDim) {
        device::Idx flat = idx[0];
        for (size_t d = 1; d < NDim; ++d)
          flat = flat * (sizesInMem[d] + 2 * nGhosts) + idx[d];

        device::Idx expectedParity = computeParity<NDim>(idx, localStarts, transForward, nGhosts);
        tdd.verify(hostParity[flat] == expectedParity);
        if (expectedParity == 0)
          ++verifiedEven;
        else
          ++verifiedOdd;
        return;
      }
      for (device::Idx i = nGhosts; i < nGhosts + sizesInMem[dim]; ++i) {
        idx[dim] = i;
        iterate(dim + 1, idx);
      }
    };

    device::IdxArray<NDim> idx{};
    iterate(0, idx);

    tdd.verify(verifiedEven == evenCount);
    tdd.verify(verifiedOdd == oddCount);
  }

  // ---- 2D transposed test ----
  template <> void CheckerboardLayoutTester<2>::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 2;
    const device::Idx nGhosts = 1;
    const device::IdxArray<NDim> gridSizes{{6, 4}};

    LayoutStruct<NDim> layout(gridSizes, nGhosts);
    layout.setLocalStarts({{2, 3}});
    layout.setTranspositionMap_memoryToGlobalSpace({{1, 0}}); // transposed

    CheckerboardLayout<NDim> even(layout, Parity::Even);
    CheckerboardLayout<NDim> odd(layout, Parity::Odd);

    device::Idx evenCount = 0, oddCount = 0;
    device::iteration::reduce(
        "CountEven2D", even, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx, device::Idx &update) { update += 1; },
        evenCount);
    device::iteration::reduce(
        "CountOdd2D", odd, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx, device::Idx &update) { update += 1; },
        oddCount);

    // After transposition: sizesInMemory[0] = localSizes[1] = 4, sizesInMemory[1] = localSizes[0] = 6
    // Total = 4 * 6 = 24
    device::Idx totalSites = 1;
    const auto sizesInMem = layout.getSizesInMemory();
    for (size_t d = 0; d < NDim; ++d)
      totalSites *= sizesInMem[d];

    say << "2D transposed: evenCount=" << evenCount << " oddCount=" << oddCount << " total=" << totalSites << "\n";

    tdd.verify(evenCount + oddCount == totalSites);
    tdd.verify(evenCount == totalSites / 2);
    tdd.verify(oddCount == totalSites / 2);
  }

  // ---- Visual test: 4x4x4 lattice with linear indices ----
  template <size_t NDim> struct CheckerboardLayoutVisualTester {
    static void Test(TDDAssertion &tdd);
  };

  template <> void CheckerboardLayoutVisualTester<3>::Test(TDDAssertion &tdd)
  {
// I DON'T UNDERSTAND WHY NVCC HAS THESE BUGS but it has bugs. constexpr on nvcc is flaky as hell, and in this test it
// just breaks totally apart
#ifndef __NVCC__
    constexpr size_t NDim = 3;
    constexpr device::Idx N = 4;
    const device::Idx nGhosts = 0;
    const device::IdxArray<NDim> gridSizes{{N, N, N}};

    LayoutStruct<NDim> layout(gridSizes, nGhosts);

    CheckerboardLayout<NDim> even(layout, Parity::Even);
    CheckerboardLayout<NDim> odd(layout, Parity::Odd);

    constexpr device::Idx totalSites = N * N * N;

    auto evenView = device::memory::NDView<device::Idx, 1>("evenVis", totalSites);
    auto oddView = device::memory::NDView<device::Idx, 1>("oddVis", totalSites);

    // Initialize to -1
    device::iteration::foreach<1>(
        "InitEvenVis", {0}, {totalSites}, DEVICE_LAMBDA(const device::IdxArray<1> &i) { evenView(i[0]) = -1; });
    device::iteration::foreach<1>(
        "InitOddVis", {0}, {totalSites}, DEVICE_LAMBDA(const device::IdxArray<1> &i) { oddView(i[0]) = -1; });

    // Fill even sites with linear index
    device::iteration::foreach (
        "FillEven", even, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx) {
          device::Idx linear = fullIdx[0] * N * N + fullIdx[1] * N + fullIdx[2];
          evenView(linear) = linear;
        });

    // Fill odd sites with linear index
    device::iteration::foreach (
        "FillOdd", odd, DEVICE_LAMBDA(const device::IdxArray<NDim> &fullIdx) {
          device::Idx linear = fullIdx[0] * N * N + fullIdx[1] * N + fullIdx[2];
          oddView(linear) = linear;
        });

    // Copy to host
    std::vector<device::Idx> hostEven(totalSites), hostOdd(totalSites);
    device::memory::copyDeviceToHost(evenView, hostEven.data());
    device::memory::copyDeviceToHost(oddView, hostOdd.data());

    // Build display as single strings so say<< prints readable blocks
    // Show spatial coordinates (x,y,z) for visited sites, --- for empty
    auto buildGrid = [&](const std::vector<device::Idx> &host, const std::string &label) {
      std::string out = "\n" + label + "\n";
      for (device::Idx i = 0; i < N; ++i) {
        out += "  Slice x=" + std::to_string(i) + ":\n";
        out += "  +---------------------------+\n";
        for (device::Idx j = 0; j < N; ++j) {
          out += "  |";
          for (device::Idx k = 0; k < N; ++k) {
            device::Idx lin = i * N * N + j * N + k;
            if (host[lin] >= 0)
              out += " (" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ")";
            else
              out += "  ---  ";
          }
          out += " |\n";
        }
        out += "  +---------------------------+\n";
      }
      return out;
    };

    say << buildGrid(hostEven, "=== EVEN sites (4x4x4) ===");
    say << buildGrid(hostOdd, "=== ODD sites (4x4x4) ===");

    // Verify: every site appears in exactly one of even or odd
    for (device::Idx lin = 0; lin < totalSites; ++lin) {
      bool inEven = (hostEven[lin] == lin);
      bool inOdd = (hostOdd[lin] == lin);
      tdd.verify(inEven != inOdd); // exactly one must be true
    }
#endif
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::CheckerboardLayoutTester<3>> test3D;
  TempLat::TDDContainer<TempLat::CheckerboardLayoutTester<2>> test2D;
  TempLat::TDDContainer<TempLat::CheckerboardLayoutVisualTester<3>> testVisual;
} // namespace
