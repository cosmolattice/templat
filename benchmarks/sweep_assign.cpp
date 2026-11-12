/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

// Cost of the sweep for assignments that read their own target at neighbouring sites (phi = f(phi)),
// against the same expression assigned to a separate field, and its dependence on the band width B
// (planes per sweep step). Used to calibrate sweepBandWidth. -DSWEEP_NGRID sets the lattice size.

#include <iostream>

#include "TempLat/session/sessionguard.h"
#include "TempLat/util/benchmark.h"

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/subtract.h"
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"

#ifndef SWEEP_NGRID
#define SWEEP_NGRID 256
#endif

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  using T = double;
  constexpr size_t nGrid = SWEEP_NGRID;
  constexpr size_t nSteps = 50;

  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, 2, false);
  toolBox->unsetVerbose();

  Field<T, NDim> phi("phi", toolBox);
  Field<T, NDim> chi("chi", toolBox);
  phi.inFourierSpace() = RandomGaussianField<T, NDim>("seed", toolBox);
  chi = phi;

  // Both variants exchange phi's ghosts every step, as the aliased one changes phi.
  auto measure = [&](Benchmark::Measurer &measurer, const std::string &tag, auto &&assign) {
    assign(); // warm-up, allocates the buffer
    for (size_t i = 0; i < nSteps; ++i)
      measurer.measure(tag, [&]() {
        phi.getMemoryManager()->setGhostsAreStale();
        assign();
        device::iteration::fence();
      });
  };
  auto r1 = [&](auto &target) { target = phi + 0.1 * LatLapl(phi); };
  auto r2 = [&](auto &target) { target = phi - 0.01 * LatLapl(LatLapl(phi)); };

  Benchmark bench([&](Benchmark::Measurer &measurer) {
    measure(measurer, "r=1 separate", [&]() { r1(chi); });
    measure(measurer, "r=1 aliased", [&]() { r1(phi); });
    measure(measurer, "r=2 separate", [&]() { r2(chi); });
    measure(measurer, "r=2 aliased", [&]() { r2(phi); });
    for (int b : {1, 2, 4, 8, 16, 32}) {
      toolBox->mSweepBuffer.setBandWidth(b);
      measure(measurer, "r=2 aliased B=" + std::to_string(b), [&]() { r2(phi); });
    }
    toolBox->mSweepBuffer.setBandWidth(0);
  });
  bench.run(1);
  bench.print();
  std::cout << "Resident device threads: " << device::iteration::concurrency() << std::endl;
}
