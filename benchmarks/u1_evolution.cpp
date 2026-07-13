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

#if defined(KICK_FUSED) || defined(KICK_SCATTER)
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/parallel/device_iteration.h"
#endif

#include <cmath>
#include <cstdio>

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  using T = double;
#ifndef U1_NGRID
#define U1_NGRID 64
#endif
#ifndef U1_NSTEPS
#define U1_NSTEPS 100
#endif
  constexpr size_t nGrid = U1_NGRID;
  constexpr size_t nGhost = 1;
  constexpr size_t nSteps = U1_NSTEPS;
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
#if defined(KICK_SCATTER)
        // Scatter kick: each plaquette computed ONCE (3 per site, vs 12 plaquette evaluations
        // for the per-direction gather below), its Im part added (+/- dt) to the four momenta
        // it touches: 2 on-site and 2 one step up. No extra field.
        //
        // The sweep is extended by one layer into the LOWER ghost shell (starts at gg-1), which
        // is what makes it correct without any reverse/additive halo: E_a(y) needs the plaquette
        // at y - b^ for every b != a, and for y on the lower face of the owned box that site is
        // exactly a lower-ghost site. The mirror-image writes that run off the upper face land in
        // E's own ghost cells, which are never read (E is only ever used site-locally, in the
        // drift and in the energy), so they are simply discarded -- and their physical content is
        // not lost: it is re-supplied by the lower-ghost plaquette on the neighbouring rank (or,
        // at one rank, by the periodic wrap of the ghost layer itself).
        //
        // Cost of correctness = recomputing the plaquettes on the ghost shell, i.e. (L+1)^3 / L^3
        // sites. Serial backend => sequential, so plain += is race-free. A threaded build would
        // additionally need a coloring (or per-thread slab ownership) for the += .
        static_assert(NDim == 3, "scatter kick prototype hardcodes NDim=3");
        auto q12 = Imag(plaq(U, 1_c, 2_c));
        auto q13 = Imag(plaq(U, 1_c, 3_c));
        auto q23 = Imag(plaq(U, 2_c, 3_c));
        GhostsHunter::apply(q12);
        GhostsHunter::apply(q13);
        GhostsHunter::apply(q23);
        auto e1 = E(1_c).getView();
        auto e2 = E(2_c).getView();
        auto e3 = E(3_c).getView();
        const auto lay = E(1_c).getLayout();
        const auto Ls = lay.getLocalSizes();
        const device::Idx gg = lay.getNGhosts();
        device::IdxArray<NDim> starts, stops;
        for (size_t d = 0; d < NDim; ++d) {
          starts[d] = gg - 1;    // one layer into the lower ghost shell
          stops[d] = gg + Ls[d]; // owned box (writes to x+1 reach the upper ghost layer)
        }
        device::iteration::foreach (
            "kick_scatter", starts, stops, DEVICE_LAMBDA(const device::IdxArray<NDim> &idx) {
              device::apply(
                  [&](auto &&...args) {
                    const device::Idx c[3] = {device::Idx(args)...};
                    const double s12 = dt * DoEval::eval(q12, args...);
                    const double s13 = dt * DoEval::eval(q13, args...);
                    const double s23 = dt * DoEval::eval(q23, args...);
                    const device::Idx u0 = c[0] + 1, u1 = c[1] + 1, u2 = c[2] + 1;
                    // plane (1,2): a=1(dim0) b=2(dim1)
                    e1(c[0], c[1], c[2]) -= s12;
                    e2(c[0], c[1], c[2]) += s12;
                    e1(c[0], u1, c[2]) += s12; // E_1(x+2^)
                    e2(u0, c[1], c[2]) -= s12; // E_2(x+1^)
                    // plane (1,3): a=1(dim0) b=3(dim2)
                    e1(c[0], c[1], c[2]) -= s13;
                    e3(c[0], c[1], c[2]) += s13;
                    e1(c[0], c[1], u2) += s13; // E_1(x+3^)
                    e3(u0, c[1], c[2]) -= s13; // E_3(x+1^)
                    // plane (2,3): a=2(dim1) b=3(dim2)
                    e2(c[0], c[1], c[2]) -= s23;
                    e3(c[0], c[1], c[2]) += s23;
                    e2(c[0], c[1], u2) += s23; // E_2(x+3^)
                    e3(c[0], u1, c[2]) -= s23; // E_3(x+2^)
                  },
                  idx);
            });
#elif defined(KICK_FUSED)
        // Sweep-fusion prototype: instead of NDim separate whole-lattice assignments
        // (each one re-streaming U from RAM), evaluate all NDim direction updates in a
        // SINGLE pass over the lattice, so U/E are streamed once per step. The RHS is
        // byte-for-byte the same expression as the per-direction version below.
        static_assert(NDim == 3, "fused kick prototype hardcodes NDim=3");
        auto kickRHS = [&](auto i) {
          return E(i) - dt * Total(j, 1, NDim,
                                   IfElse(i != j, Imag(plaq(U, i, j)) - shift(Imag(plaq(U, i, j)), -j), ZeroType()));
        };
        auto r1 = kickRHS(1_c);
        auto r2 = kickRHS(2_c);
        auto r3 = kickRHS(3_c);
        // Refresh U's ghosts exactly as the per-assignment path would (assign() calls
        // GhostsHunter::apply on its RHS); the hand-written foreach below bypasses it.
        GhostsHunter::apply(r1);
        GhostsHunter::apply(r2);
        GhostsHunter::apply(r3);
        auto v1 = E(1_c).getView();
        auto v2 = E(2_c).getView();
        auto v3 = E(3_c).getView();
        device::iteration::foreach ("kick_fused", E(1_c).getLayout(),
                                    DEVICE_LAMBDA(const device::IdxArray<NDim> &idx) {
                                      device::apply(
                                          [&](auto &&...args) {
#ifdef KICK_FUSED_TMP
                                            // All reads first, then all writes: removes the
                                            // store-between-loads aliasing barrier so the
                                            // compiler MAY CSE the U-link loads shared across
                                            // directions.
                                            const auto t1 = DoEval::eval(r1, args...);
                                            const auto t2 = DoEval::eval(r2, args...);
                                            const auto t3 = DoEval::eval(r3, args...);
                                            v1(args...) = t1;
                                            v2(args...) = t2;
                                            v3(args...) = t3;
#else
                                            v1(args...) = DoEval::eval(r1, args...);
                                            v2(args...) = DoEval::eval(r2, args...);
                                            v3(args...) = DoEval::eval(r3, args...);
#endif
                                          },
                                          idx);
                                    });
        // E is only ever read site-locally (drift, energy), so its ghosts are never
        // needed; no setGhostsAreStale() required here.
#else
        for_in_range<1, NDim + 1>([&](auto i) {
          E(i) = E(i) - dt * Total(j, 1, NDim,
                                   IfElse(i != j, Imag(plaq(U, i, j)) - shift(Imag(plaq(U, i, j)), -j), ZeroType()));
        });
#endif
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
