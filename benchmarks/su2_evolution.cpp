/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

/* Pure-gauge SU(2) Yang-Mills Hamiltonian evolution in 3+1 dimensions (temporal gauge).
 * Same physics and same deterministic initial condition as bench-hila-templat's bench_su2.cpp,
 * so the energy drift below is directly comparable and acts as an exact correctness gate.
 *
 * Kick-drift leapfrog:
 *   kick : Pi_i -= dt * sum_{j != i} [ plaq(U,i,j) - plaqBack(U,i,j) ]   (projected onto the algebra)
 *   drift: U_i   = exp(dt * Pi_i) * U_i                                  (exact on the group)
 *
 * THE SCATTER KICK (-DKICK_SCATTER)
 * ---------------------------------
 * The gather kick above evaluates, per site, 3 directions x 2 transverse j's x (plaq + plaqBack)
 * = 12 four-link products. But the four plaquettes touching a given plane are all the SAME loop:
 *
 *   (1) P_ba(x) = P_ab(x)^dagger                                  [reversed traversal]
 *   (2) plaqBack_ab(x) = U_b(x-b)^dagger . P_ab(x-b) . U_b(x-b)   [see plaquetteback.h]
 *
 * (2) says the backward plaquette is a CONJUGATION of a forward one, not a new loop. For SU(2)
 * conjugation acts on the algebra 3-vector as an SO(3) rotation (the adjoint action) and leaves
 * the identity component alone -- so it costs a vector rotation, not a matrix product.
 *
 * Hence: visit each site once, compute the 3 unique plaquettes P_12, P_13, P_23 (3 four-link
 * products instead of 12), and scatter each one's algebra part p to the four momenta of its four
 * links. With A = U_a(x), D = U_b(x):
 *
 *   Pi_a(x)   -= dt * p              Pi_a(x+b) += dt * Ad(D^dagger) p
 *   Pi_b(x)   += dt * p              Pi_b(x+a) -= dt * Ad(A^dagger) p
 *
 * which reduces to exactly the U(1) scatter (Ad -> identity, p -> Im P) as it must.
 *
 * The sweep runs one layer into the LOWER ghost shell so that every owned momentum receives its
 * x-b neighbour term locally; the mirror-image writes that run off the upper face land in Pi's own
 * ghost cells, which are never read (Pi is only used site-locally, in the drift and the energy).
 * That is what makes the scatter correct under MPI with NO reverse/additive halo and no extra
 * communication. Serial backend => sequential => the += needs no atomics; a threaded build would
 * need a coloring or per-thread slab ownership.
 *
 * -DKICK_SCATTER_QMUL evaluates Ad() with two explicit SU(2) quaternion products instead of the
 * Rodrigues rotation -- slower, but an independent check of the rotation's sign convention.
 */

#include "TempLat/session/sessionguard.h"
#include "TempLat/util/benchmark.h"

#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/algebra/su2algebra/su2liealgebrafield.h"
#include "TempLat/lattice/algebra/su2algebra/su2expmap.h"
#include "TempLat/lattice/algebra/su2algebra/su2multiply.h"
#include "TempLat/lattice/algebra/su2algebra/su2sum.h"
#include "TempLat/lattice/algebra/su2algebra/su2subtract.h"
#include "TempLat/lattice/algebra/su2algebra/scalarsu2multiplication.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquette.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquetteback.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/staticif.h"

#if defined(KICK_SCATTER)
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/util/constexpr_for.h"
#endif

#include <cmath>
#include <cstdio>

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  using T = double;
#ifndef SU2_NGRID
#define SU2_NGRID 64
#endif
#ifndef SU2_NSTEPS
#define SU2_NSTEPS 100
#endif
  constexpr size_t nGrid = SU2_NGRID;
  constexpr size_t nGhost = 1;
  constexpr size_t nSteps = SU2_NSTEPS;
  constexpr T dt = 0.01;
  constexpr T amp = 0.5;

  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost, false);
  toolBox->unsetVerbose();

  FieldCollection<SU2Field<T, NDim>, NDim, false, 1> U("U", toolBox);
  FieldCollection<SU2LieAlgebraField<T, NDim>, NDim, false, 1> Pi("Pi", toolBox);

  // Deterministic plane-wave links U_i = exp(i amp sin(2 pi x_t / N) sigma_i), Pi = 0.
  // Each angle varies transverse to its link and uses a different Pauli matrix, so the
  // configuration is genuinely non-abelian and not pure gauge.
  SpatialCoordinate<NDim> x(toolBox);
  const T w = 2.0 * M_PI / nGrid;

  U(1_c)(0_c) = cos(amp * sin(w * getVectorComponent(x, Tag<1>())));
  U(1_c)(1_c) = sin(amp * sin(w * getVectorComponent(x, Tag<1>())));
  U(1_c)(2_c) = T(0);
  U(1_c)(3_c) = T(0);
  U(2_c)(0_c) = cos(amp * sin(w * getVectorComponent(x, Tag<2>())));
  U(2_c)(1_c) = T(0);
  U(2_c)(2_c) = sin(amp * sin(w * getVectorComponent(x, Tag<2>())));
  U(2_c)(3_c) = T(0);
  U(3_c)(0_c) = cos(amp * sin(w * getVectorComponent(x, Tag<0>())));
  U(3_c)(1_c) = T(0);
  U(3_c)(2_c) = T(0);
  U(3_c)(3_c) = sin(amp * sin(w * getVectorComponent(x, Tag<0>())));

  for_in_range<1, NDim + 1>([&](auto i) {
    Pi(i)(1_c) = T(0);
    Pi(i)(2_c) = T(0);
    Pi(i)(3_c) = T(0);
  });

  // H = 2 sum_i <Pi_i . Pi_i> + 2 sum_{i != j} <1 - c0(P_ij)>, per site.
  auto energy = [&]() {
    T elec = 2.0 * Total(i, 1, NDim, average(pow<2>(Pi(i)(1_c)) + pow<2>(Pi(i)(2_c)) + pow<2>(Pi(i)(3_c))));
    T mag = 2.0 * Total(i, 1, NDim, Total(j, 1, NDim, IfElse(i != j, average(1.0 - plaq(U, i, j).SU2Get(0_c)), 0.0)));
    return std::make_pair(elec, mag);
  };

  const auto [elec0, mag0] = energy();

  Benchmark bench([&](Benchmark::Measurer &measurer) {
    for (size_t step = 0; step < nSteps; ++step) {
      measurer.measure("Pi kick", [&]() {
#if defined(KICK_SCATTER)
        static_assert(NDim == 3, "scatter kick prototype hardcodes NDim=3");

        auto p12 = plaq(U, 1_c, 2_c);
        auto p13 = plaq(U, 1_c, 3_c);
        auto p23 = plaq(U, 2_c, 3_c);
        GhostsHunter::apply(p12);
        GhostsHunter::apply(p13);
        GhostsHunter::apply(p23);

        // Link views: u[a][k] = component k of U_a. Momentum views: pi[a][k], k = 0..2 <-> algebra 1..3.
        const auto u1 = U(1_c), u2 = U(2_c), u3 = U(3_c);
        const device::memory::NDViewUnmanaged<T, NDim> uv[3][4] = {
            {u1(0_c).getView(), u1(1_c).getView(), u1(2_c).getView(), u1(3_c).getView()},
            {u2(0_c).getView(), u2(1_c).getView(), u2(2_c).getView(), u2(3_c).getView()},
            {u3(0_c).getView(), u3(1_c).getView(), u3(2_c).getView(), u3(3_c).getView()}};
        const auto q1 = Pi(1_c), q2 = Pi(2_c), q3 = Pi(3_c);
        const device::memory::NDViewUnmanaged<T, NDim> pv[3][3] = {
            {q1(1_c).getView(), q1(2_c).getView(), q1(3_c).getView()},
            {q2(1_c).getView(), q2(2_c).getView(), q2(3_c).getView()},
            {q3(1_c).getView(), q3(2_c).getView(), q3(3_c).getView()}};

        const auto lay = q1(1_c).getLayout();
        const auto Ls = lay.getLocalSizes();
        const device::Idx gg = lay.getNGhosts();
        device::IdxArray<NDim> starts, stops;
        for (size_t d = 0; d < NDim; ++d) {
          starts[d] = gg - 1;    // one layer into the lower ghost shell
          stops[d] = gg + Ls[d]; // writes to x+1 reach the upper ghost layer
        }

        device::iteration::foreach (
            "su2_kick_scatter", starts, stops, DEVICE_LAMBDA(const device::IdxArray<NDim> &idx) {
              device::apply(
                  [&](auto &&...args) {
                    const device::Idx i0 = idx[0], i1 = idx[1], i2 = idx[2];
                    const device::Idx j0 = i0 + 1, j1 = i1 + 1, j2 = i2 + 1;

                    // Rotate the algebra 3-vector p by Ad(G^dagger), G a unit quaternion (g0, g):
                    //   Ad(G^dagger) p = (g0^2 - |g|^2) p + 2 (g.p) g + 2 g0 (g x p)
                    // (derived from the quaternion product convention in su2multiply.h).
                    // KICK_SCATTER_QMUL instead does the two explicit quaternion products
                    // Q = (G^dagger . P) . G with vec(r,t) = r0 t + t0 r - r x t, as a sign check.
                    auto adDagger = [](const T g[4], const T p[3], T out[3]) {
#if defined(KICK_SCATTER_QMUL)
                      const T gd[4] = {g[0], -g[1], -g[2], -g[3]};
                      const T M0 = -(gd[1] * p[0] + gd[2] * p[1] + gd[3] * p[2]);
                      const T Mv[3] = {gd[0] * p[0] - (gd[2] * p[2] - gd[3] * p[1]),
                                       gd[0] * p[1] - (gd[3] * p[0] - gd[1] * p[2]),
                                       gd[0] * p[2] - (gd[1] * p[1] - gd[2] * p[0])};
                      out[0] = M0 * g[1] + g[0] * Mv[0] - (Mv[1] * g[3] - Mv[2] * g[2]);
                      out[1] = M0 * g[2] + g[0] * Mv[1] - (Mv[2] * g[1] - Mv[0] * g[3]);
                      out[2] = M0 * g[3] + g[0] * Mv[2] - (Mv[0] * g[2] - Mv[1] * g[1]);
#else
                      const T g0 = g[0], g1 = g[1], g2 = g[2], g3 = g[3];
                      const T s = g0 * g0 - (g1 * g1 + g2 * g2 + g3 * g3);
                      const T d2 = T(2) * (g1 * p[0] + g2 * p[1] + g3 * p[2]);
                      const T t0 = T(2) * g0;
                      out[0] = s * p[0] + d2 * g1 + t0 * (g2 * p[2] - g3 * p[1]);
                      out[1] = s * p[1] + d2 * g2 + t0 * (g3 * p[0] - g1 * p[2]);
                      out[2] = s * p[2] + d2 * g3 + t0 * (g1 * p[1] - g2 * p[0]);
#endif
                    };

                    // The three link quaternions at x, read ONCE (each is used by two planes).
                    const T L1[4] = {uv[0][0](i0, i1, i2), uv[0][1](i0, i1, i2), uv[0][2](i0, i1, i2),
                                     uv[0][3](i0, i1, i2)};
                    const T L2[4] = {uv[1][0](i0, i1, i2), uv[1][1](i0, i1, i2), uv[1][2](i0, i1, i2),
                                     uv[1][3](i0, i1, i2)};
                    const T L3[4] = {uv[2][0](i0, i1, i2), uv[2][1](i0, i1, i2), uv[2][2](i0, i1, i2),
                                     uv[2][3](i0, i1, i2)};

                    // One plane (a, b), fully unrolled: a/b are template params so every view index
                    // and every neighbour offset folds at compile time.
                    auto plane = [&]<int a, int b>(const device::array<T, 4> &P, const T A[4], const T D[4],
                                                   device::Idx ba0, device::Idx ba1, device::Idx ba2, device::Idx ab0,
                                                   device::Idx ab1, device::Idx ab2) {
                      const T p[3] = {dt * P[1], dt * P[2], dt * P[3]};
                      T rD[3], rA[3];
                      adDagger(D, p, rD); // for Pi_a(x + b^)
                      adDagger(A, p, rA); // for Pi_b(x + a^)
                      constexpr_for<0, 3>([&](auto m) {
                        pv[a][m](i0, i1, i2) -= p[m];     // Pi_a(x)     -= dt p
                        pv[b][m](i0, i1, i2) += p[m];     // Pi_b(x)     += dt p
                        pv[a][m](ba0, ba1, ba2) += rD[m]; // Pi_a(x + b^) += dt Ad(D^dag) p
                        pv[b][m](ab0, ab1, ab2) -= rA[m]; // Pi_b(x + a^) -= dt Ad(A^dag) p
                      });
                    };

                    // plane (1,2): a = dim0, b = dim1
                    plane.template operator()<0, 1>(DoEval::eval(p12, args...), L1, L2, i0, j1, i2, j0, i1, i2);
                    // plane (1,3): a = dim0, b = dim2
                    plane.template operator()<0, 2>(DoEval::eval(p13, args...), L1, L3, i0, i1, j2, j0, i1, i2);
                    // plane (2,3): a = dim1, b = dim2
                    plane.template operator()<1, 2>(DoEval::eval(p23, args...), L2, L3, i0, i1, j2, i0, j1, i2);
                  },
                  idx);
            });
#elif defined(KICK_FUSED)
        // The natural way to write it -- and 12-20% slower. The whole Total over j is one expression,
        // so all 12 four-link products are live at once; vectorized, that is far more than the 32 zmm
        // registers hold, and GCC spills (3276 stack-touching vector moves in the loop, against 1500
        // for the split form below and 83 for a hand-written kernel). Kept so the comparison stays
        // reproducible: benchmarks/PERFORMANCE.md section 2.
        for_in_range<1, NDim + 1>([&](auto i) {
          Pi(i) = Pi(i) - dt * Total(j, 1, NDim,
                                     IfElse(i != j, plaq(U, i, j) - plaqBack(U, i, j), ZeroType()));
        });
#else
        // Same arithmetic, accumulated one transverse direction at a time. This costs one extra Pi
        // read+write per j and buys a live set that fits in the register file -- a net 12-20%. Do NOT
        // split further (one assignment per term): the extra Pi traffic then outweighs the spills it
        // saves. Energy drift is bit-identical to the fused form, as it must be.
        for_in_range<1, NDim + 1>([&](auto i) {
          for_in_range<1, NDim + 1>([&](auto j) {
            if constexpr (decltype(i)::value != decltype(j)::value) {
              Pi(i) = Pi(i) - dt * (plaq(U, i, j) - plaqBack(U, i, j));
            }
          });
        });
#endif
        device::iteration::fence();
      });
      measurer.measure("U drift", [&]() {
        for_in_range<1, NDim + 1>([&](auto i) { U(i) = exp(dt * Pi(i)) * U(i); });
        device::iteration::fence();
      });
    }
  });
  bench.run(1);
  bench.print();
  bench.log("su2_evolution");

  const auto [elec1, mag1] = energy();
  std::printf("initial: electric = %.8f  magnetic = %.8f  total = %.8f\n", elec0, mag0, elec0 + mag0);
  std::printf("final:   electric = %.8f  magnetic = %.8f  total = %.8f\n", elec1, mag1, elec1 + mag1);
  std::printf("plaq = %.10g\n", (double)(mag1 / 12.0));
  std::printf("relative energy drift after %zu steps: %.3e\n", nSteps, std::abs((elec1 + mag1) / (elec0 + mag0) - 1.0));
}
