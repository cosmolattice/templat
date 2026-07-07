/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

/* Compact U(1) pure-gauge Hamiltonian evolution in 3+1 dimensions (temporal gauge).
 *
 * Degrees of freedom:
 *   - links     U_i(x) = e^{i theta_i(x)}  -> one U1Field (= ComplexField) per direction
 *   - electric  E_i(x)                     -> one real Field per direction
 *
 * Hamiltonian (lattice units, dx = 1, coupling absorbed):
 *   H = 1/2 sum_{x,i} E_i(x)^2 + sum_{x,i<j} [ 1 - Re P_ij(x) ]
 * with the plaquette P_ij(x) = U_i(x) U_j(x+i) U_i*(x+j) U_j*(x).
 *
 * Equations of motion, integrated with a kick-drift leapfrog:
 *   dE_i/dt = -sum_{j != i} Im[ P_ij(x) - P_ij(x - j) ]
 *   dU_i/dt = i E_i U_i   ->  multiplicative update, exact on the group:
 *                             U_i <- complexPhase(dt * E_i) * U_i
 *
 * The total energy is printed before and after the evolution; it is conserved up to the
 * O(dt^2) leapfrog error while energy oscillates between the electric and magnetic sectors.
 */

#include "TempLat/session/sessionguard.h"
#include "TempLat/util/benchmark.h"

#include "TempLat/lattice/algebra/gaugealgebra/u1field.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquette.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/staticif.h"

#include <cmath>
#include <cstdio>

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  using T = double;
  constexpr size_t nGrid = 64;
  constexpr size_t nGhost = 1;
  constexpr size_t nSteps = 100;
  constexpr T dt = 0.02;

  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost, false);
  toolBox->unsetVerbose();

  // Gauge links U_1..U_NDim and electric fields E_1..E_NDim, indexed by direction tags 1..NDim.
  FieldCollection<U1Field<T, NDim>, NDim, false, 1> U("U", toolBox);
  FieldCollection<Field<T, NDim>, NDim, false, 1> E("E", toolBox);

  // Initial conditions: plane-wave angles theta_i(x) = 0.5 sin(2 pi x_{i mod NDim} / N), E = 0.
  // Each theta_i varies along a direction other than i, so the initial magnetic energy is nonzero
  // (theta_i varying only along direction i would be a pure-gauge configuration with F = 0).
  SpatialCoordinate x(toolBox);
  for_in_range<1, NDim + 1>([&](auto i) {
    U(i) = complexPhase(0.5 * sin((2.0 * M_PI / nGrid) * getVectorComponent(x, Tag<i % NDim>())));
    E(i) = T(0);
  });

  // H = 1/2 sum_i <E_i^2> + sum_{i<j} <1 - Re P_ij>, per site.
  auto energy = [&]() {
    T elec = 0.5 * Total(i, 1, NDim, average(pow<2>(E(i))));
    T mag = 0.5 * Total(i, 1, NDim, Total(j, 1, NDim, IfElse(i != j, average(1.0 - Real(plaq(U, i, j))), 0.0)));
    return std::make_pair(elec, mag);
  };

  const auto [elec0, mag0] = energy();

  Benchmark bench([&](Benchmark::Measurer &measurer) {
    for (size_t step = 0; step < nSteps; ++step) {
      // Kick: dE_i/dt = -sum_{j != i} Im[ P_ij(x) - P_ij(x - j) ]
      measurer.measure("E kick", [&]() {
        for_in_range<1, NDim + 1>([&](auto i) {
          E(i) = E(i) - dt * Total(j, 1, NDim,
                                   IfElse(i != j, Imag(plaq(U, i, j)) - shift(Imag(plaq(U, i, j)), -j), ZeroType()));
        });
        device::iteration::fence();
      });
      // Drift: U_i <- e^{i dt E_i} U_i. Multiplicative update, stays exactly on the group,
      // so no unitarization is ever needed. Aliasing U on both sides is safe: the right
      // hand side is site-local (the kick/drift split keeps shifts of U out of U's own update).
      measurer.measure("U drift", [&]() {
        for_in_range<1, NDim + 1>([&](auto i) { U(i) = complexPhase(dt * E(i)) * U(i); });
        device::iteration::fence();
      });
    }
  });
  bench.run(1);
  bench.print();
  bench.log("u1_evolution");

  const auto [elec1, mag1] = energy();
  std::printf("initial: electric = %.8f  magnetic = %.8f  total = %.8f\n", elec0, mag0, elec0 + mag0);
  std::printf("final:   electric = %.8f  magnetic = %.8f  total = %.8f\n", elec1, mag1, elec1 + mag1);
  std::printf("relative energy drift after %zu steps: %.3e\n", nSteps, std::abs((elec1 + mag1) / (elec0 + mag0) - 1.0));
}
