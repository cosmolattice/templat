/* Diagnostic benchmark for Kokkos backend threading performance.
 *
 * Purpose: isolate the root causes of the OpenMP vs Threads performance gap
 * by benchmarking progressively more complex kernels:
 *
 *  1. fill_const   — pure store (memory bandwidth / thread overhead baseline)
 *  2. copy         — field-to-field copy (read + write)
 *  3. add          — simple expression without stencil (no ghost accesses)
 *  4. init_gauss   — Gaussian random field assignment in Fourier space
 *                    (reproduces the initialize_field anomaly from phi4)
 *  5. raw_tile_N   — direct Kokkos MDRangePolicy sweep over tile sizes
 *                    {1,8,16,32,64,128,256}, bypassing TempLat field
 *                    machinery, to isolate tiling granularity effects.
 *
 * Build both OpenMP and Threads builds and compare per-operation timings.
 * For the raw tile sweep, look for the tile size where throughput peaks —
 * this is the optimal value for TL_TILE_SIZE_INNER on this machine.
 *
 * NUMA diagnostic: if initialize_field is anomalously slow on Threads,
 * rerun with `numactl --interleave=all ./bench-threading_perf` to check
 * whether NUMA first-touch is the culprit.
 */

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/session/sessionguard.h"
#include "TempLat/util/benchmark.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/parallel/devices/kokkos/kokkos_memory.h"

#include <Kokkos_Core.hpp>
#include <string>

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  using T = double;
  constexpr size_t nGrid = 256;
  constexpr size_t nGhost = 1;

  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost, false);
  toolBox->unsetVerbose();

  Field<T, NDim> phi("phi", toolBox);
  Field<T, NDim> phi2("phi2", toolBox);

  // Touch all config-space allocations before benchmarking
  phi.getMemoryManager()->confirmConfigSpace();
  phi2.getMemoryManager()->confirmConfigSpace();
  phi.getMemoryManager()->confirmFourierSpace();
  phi2.getMemoryManager()->confirmFourierSpace();
  phi.getMemoryManager()->confirmConfigSpace();
  phi2.getMemoryManager()->confirmConfigSpace();

  Benchmark bench([&](Benchmark::Measurer &measurer) {

    // ------------------------------------------------------------------
    // Kernel 1: fill with a constant — pure write, no compute
    // ------------------------------------------------------------------
    measurer.measure("fill_const", [&]() {
      phi = T(1.0);
      device::iteration::fence();
    });

    // ------------------------------------------------------------------
    // Kernel 2: field copy — read + write
    // ------------------------------------------------------------------
    measurer.measure("copy", [&]() {
      phi2 = phi;
      device::iteration::fence();
    });

    // ------------------------------------------------------------------
    // Kernel 3: simple expression, no stencil, no ghost access
    // ------------------------------------------------------------------
    measurer.measure("scale_add", [&]() {
      phi2 = T(2.0) * phi;
      device::iteration::fence();
    });

    // ------------------------------------------------------------------
    // Kernel 4: Gaussian random field in Fourier space
    // Reproduces the phi4 "initialize_field" operation.
    // ------------------------------------------------------------------
    measurer.measure("init_gauss", [&]() {
      phi.inFourierSpace() = RandomGaussianField<T, NDim>("rng", toolBox);
      device::iteration::fence();
    });

    // ------------------------------------------------------------------
    // Kernel 5: direct MDRangePolicy tile sweep
    // Bypasses TempLat field machinery. Allocates a raw Kokkos View and
    // runs a pure-store kernel with each tile size. The tile size where
    // throughput peaks is the optimal TL_TILE_SIZE_INNER for this machine.
    //
    // Note: the View is allocated once outside the measured loop. The loop
    // over tile sizes is *inside* the measurer so each tile gets its own
    // timing entry. A fence is placed before and after (the Measurer class
    // already fences on entry).
    // ------------------------------------------------------------------
    {
      constexpr size_t viewSize = nGrid + 2 * nGhost;
      device_kokkos::memory::NDView<double, NDim>
          rawView("threading_perf_raw", viewSize, viewSize, viewSize);

      // Warm up the view with a plain fill so pages are touched before timing
      Kokkos::deep_copy(rawView, T(0));
      device::iteration::fence();

      constexpr int64_t lo = static_cast<int64_t>(nGhost);
      constexpr int64_t hi = static_cast<int64_t>(nGhost + nGrid);

      for (int64_t tile_inner : {int64_t{1}, int64_t{8}, int64_t{16},
                                 int64_t{32}, int64_t{64}, int64_t{128},
                                 int64_t{256}}) {
        const std::string tag = "raw_tile_" + std::to_string(tile_inner);
        measurer.measure(tag, [&]() {
          Kokkos::Array<int64_t, 3> start = {{lo, lo, lo}};
          Kokkos::Array<int64_t, 3> stop  = {{hi, hi, hi}};
          Kokkos::Array<int64_t, 3> tiles = {{1, 1, tile_inner}};
          Kokkos::parallel_for(
              tag,
              Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace,
                                    Kokkos::Rank<3>>(start, stop, tiles),
              KOKKOS_LAMBDA(int64_t i, int64_t j, int64_t k) {
                rawView(i, j, k) = T(1);
              });
          device::iteration::fence();
        });
      }
    }
  });

  bench.run(1);
  bench.print();
  bench.log("threading_perf");
}
