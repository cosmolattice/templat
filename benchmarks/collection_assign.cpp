/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

// Assigning a FieldCollection of N fields in one kernel (fc = expr) against one kernel per component, for an
// O(N) drift and kick. -DCOLL_NGRID sets the lattice size.

#include "TempLat/session/sessionguard.h"
#include "TempLat/util/benchmark.h"

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/subtract.h"
#include "TempLat/lattice/algebra/listoperators/listoperators.h"
#include "TempLat/lattice/algebra/listoperators/norm.h"
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/util/rangeiteration/for_in_range.h"

#ifndef COLL_NGRID
#define COLL_NGRID 256
#endif

namespace TempLat
{
  template <typename C, typename R> void perComponent(C &target, R &&r)
  {
    using nakedR = std::remove_reference_t<R>;
    for_in_range<0, nakedR::size>([&](auto i) { target.getComp(i) = nakedR::Getter::get(r, i); });
  }

  template <int N>
  void run(Benchmark::Measurer &measurer, const device::memory::host_ptr<MemoryToolBox<3>> &toolBox, size_t nSteps)
  {
    using T = double;
    FieldCollection<Field<T, 3>, N> phi("phi", toolBox), pi("pi", toolBox);
    for_in_range<0, N>([&](auto i) {
      phi(i).inFourierSpace() = RandomGaussianField<T, 3>("seed" + std::to_string(int(i)), toolBox);
      pi(i) = 0.1 * phi(i);
    });
    const T dt = 0.01, lambda = 0.1;

    auto measure = [&](const std::string &tag, auto &&assign) {
      assign(); // warm-up
      for (size_t s = 0; s < nSteps; ++s)
        measurer.measure("N=" + std::to_string(N) + " " + tag, [&]() {
          assign();
          device::iteration::fence();
        });
    };
    measure("drift per-component", [&]() { perComponent(phi, phi + dt * pi); });
    measure("drift one kernel   ", [&]() { phi = phi + dt * pi; });
    measure("kick  per-component", [&]() { perComponent(pi, pi + dt * (LatLapl(phi) - lambda * norm2(phi) * phi)); });
    measure("kick  one kernel   ", [&]() { pi = pi + dt * (LatLapl(phi) - lambda * norm2(phi) * phi); });
  }
} // namespace TempLat

int main(int argc, char **argv)
{
  using namespace TempLat;
  SessionGuard guard(argc, argv, false);
  constexpr size_t nGrid = COLL_NGRID;
  constexpr size_t nSteps = 20;

  auto toolBox = MemoryToolBox<3>::makeShared(nGrid, 1, false);
  toolBox->unsetVerbose();

  Benchmark bench([&](Benchmark::Measurer &measurer) {
    run<2>(measurer, toolBox, nSteps);
    run<4>(measurer, toolBox, nSteps);
    run<8>(measurer, toolBox, nSteps);
  });
  bench.run(1);
  bench.print();
}
