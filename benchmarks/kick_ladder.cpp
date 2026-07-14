/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

/* A ladder of SU(2) kicks, from TempLat's own kernel to a hand-written ceiling, to localise where
 * the ~2.3x of the scalar-codegen gap actually lives. All variants do the SAME arithmetic on the
 * SAME 12-array-SoA layout; they differ only in how the site loop and the leaf loads are expressed.
 *
 *   -DKICK_VARIANT=0  TempLat as-is: Pi(i) = Pi(i) - dt * Total(j, plaq - plaqBack)
 *   -DKICK_VARIANT=1  same expression tree, but the site loop is a hand-written triple loop here
 *                     calling DoEval::eval(e, x, y, z)   -> isolates Kokkos' MDRange dispatch
 *   -DKICK_VARIANT=2  as 1, but the leaves load through raw T* __restrict__ instead of
 *                     Kokkos::View::operator()  -> isolates the View accessor. MEASURED AND REFUTED:
 *                     identical to the View to within noise. The code is gone; the rung #errors.
 *   -DKICK_VARIANT=3  the ceiling: plain restrict-qualified triple loop, no expression templates,
 *                     no Kokkos. Performance probe only (own arrays, not bit-compared).
 *   -DKICK_VARIANT=4  the fused kick split into one assignment per transverse direction j. The fused
 *                     expression's live set does not fit in the register file, so it SPILLS -- this
 *                     is what the remaining gap to the ceiling is made of. 12-20% faster than 0.
 *   -DKICK_VARIANT=5  split further, one assignment per term. Over-splitting: the extra Pi traffic
 *                     costs more than the spills it saves. Slower than 4.
 *
 * Variants 0, 1, 4, 5 share su2_evolution.cpp's deterministic plane-wave IC, so plaq and the relative
 * energy drift are an exact correctness gate and must be identical across all of them (they are).
 */

#include "TempLat/session/sessionguard.h"

#include "TempLat/lattice/algebra/gaugealgebra/plaquette.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquetteback.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/su2algebra/scalarsu2multiplication.h"
#include "TempLat/lattice/algebra/su2algebra/su2expmap.h"
#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/algebra/su2algebra/su2liealgebrafield.h"
#include "TempLat/lattice/algebra/su2algebra/su2multiply.h"
#include "TempLat/lattice/algebra/su2algebra/su2subtract.h"
#include "TempLat/lattice/algebra/su2algebra/su2sum.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/staticif.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef KICK_VARIANT
#define KICK_VARIANT 0
#endif
#ifndef SU2_NGRID
#define SU2_NGRID 64
#endif
#ifndef SU2_NSTEPS
#define SU2_NSTEPS 100
#endif

static constexpr size_t nGrid = SU2_NGRID;
static constexpr size_t nSteps = SU2_NSTEPS;
static constexpr double dt = 0.01;

// ---------------------------------------------------------------------------------------------
// Variant 3: the ceiling. Same arithmetic, same SoA layout, but a plain restrict-qualified loop.
// ---------------------------------------------------------------------------------------------
#if KICK_VARIANT == 3

#include <type_traits>
#include <utility>

static constexpr int N = (int)nGrid;
static constexpr int M = N + 2; // ghost-padded extent, nGhost = 1
static constexpr long SZ = (long)M * M * M;
static constexpr long S_[3] = {(long)M * M, (long)M, 1L};

using R = double *__restrict__;

static inline void qmul(const double a0, const double a1, const double a2, const double a3, const double b0,
                        const double b1, const double b2, const double b3, double &c0, double &c1, double &c2,
                        double &c3)
{
  c0 = a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3;
  c1 = a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2;
  c2 = a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1;
  c3 = a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0;
}

int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  static double *U[3][4], *P[3][3];
  for (int d = 0; d < 3; ++d) {
    for (int c = 0; c < 4; ++c) {
      U[d][c] = (double *)aligned_alloc(64, SZ * sizeof(double));
      for (long i = 0; i < SZ; ++i)
        U[d][c][i] = (c == 0) ? 1.0 : 0.001 * (d + c + (i % 7));
    }
    for (int c = 0; c < 3; ++c) {
      P[d][c] = (double *)aligned_alloc(64, SZ * sizeof(double));
      std::memset(P[d][c], 0, SZ * sizeof(double));
    }
  }

  constexpr long sx = (long)M * M, sy = M;

  // i and j are compile-time in TempLat (Tag<>), so the z-loop body must be straight-line.
  auto kick_dir = [&]<int i>(std::integral_constant<int, i>) {
    R pi0 = P[i][0], pi1 = P[i][1], pi2 = P[i][2];
    for (int x = 1; x <= N; ++x)
      for (int y = 1; y <= N; ++y) {
        const long base = (long)x * sx + (long)y * sy;
#pragma GCC ivdep
        for (int z = 1; z <= N; ++z) {
          const long q = base + z;
          double f1 = 0, f2 = 0, f3 = 0;
          [&]<int... js>(std::integer_sequence<int, js...>) {
            (
                [&] {
                  constexpr int j = js;
                  if constexpr (j != i) {
                    constexpr long si = S_[i], sj = S_[j];
                    double a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3;
                    // forward plaquette  P_ij(q) = U_i(q) U_j(q+i) U_i(q+j)^t U_j(q)^t
                    qmul(U[i][0][q], U[i][1][q], U[i][2][q], U[i][3][q], U[j][0][q + si], U[j][1][q + si],
                         U[j][2][q + si], U[j][3][q + si], a0, a1, a2, a3);
                    qmul(U[i][0][q + sj], -U[i][1][q + sj], -U[i][2][q + sj], -U[i][3][q + sj], U[j][0][q], -U[j][1][q],
                         -U[j][2][q], -U[j][3][q], b0, b1, b2, b3);
                    qmul(a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3);
                    f1 += c1;
                    f2 += c2;
                    f3 += c3;
                    // backward plaquette at q - j
                    const long r = q - sj;
                    qmul(U[i][0][r], U[i][1][r], U[i][2][r], U[i][3][r], U[j][0][r + si], U[j][1][r + si],
                         U[j][2][r + si], U[j][3][r + si], a0, a1, a2, a3);
                    qmul(U[i][0][r + sj], -U[i][1][r + sj], -U[i][2][r + sj], -U[i][3][r + sj], U[j][0][r], -U[j][1][r],
                         -U[j][2][r], -U[j][3][r], b0, b1, b2, b3);
                    qmul(a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3);
                    f1 -= c1;
                    f2 -= c2;
                    f3 -= c3;
                  }
                }(),
                ...);
          }(std::make_integer_sequence<int, 3>{});
          pi0[q] -= dt * f1;
          pi1[q] -= dt * f2;
          pi2[q] -= dt * f3;
        }
      }
  };

  auto t0 = std::chrono::steady_clock::now();
  for (size_t step = 0; step < nSteps; ++step) {
    kick_dir(std::integral_constant<int, 0>{});
    kick_dir(std::integral_constant<int, 1>{});
    kick_dir(std::integral_constant<int, 2>{});
  }
  auto t1 = std::chrono::steady_clock::now();

  const double s = std::chrono::duration<double>(t1 - t0).count() / nSteps;
  double sink = 0;
  for (int d = 0; d < 3; ++d)
    sink += P[d][0][SZ / 2] + P[d][1][SZ / 3];
  std::printf("variant  = 3 (hand-written ceiling)\n");
  std::printf("N        = %d\n", N);
  std::printf("kick_ns  = %.1f\n", 1e9 * s / ((double)N * N * N));
  std::printf("drift_ns = 0.0\n");
  std::printf("plaq     = n/a (sink %.6g)\n", sink);
  std::printf("edrift   = n/a\n");
  return 0;
}

#else // ---------------------------------------------------------- variants 0, 1, 2 (TempLat)

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  using T = double;
  constexpr size_t nGhost = 1;
  constexpr T amp = 0.5;

  // A cubic box CANNOT see the cache-blocking regression that the cluster hit: its contiguous rows are
  // short, so sweeping a whole row per (x,y) keeps the stencil's neighbour rows in cache no matter
  // what the dispatch does. The 2D-pencil decomposition MPI actually uses leaves the contiguous
  // dimension UNDIVIDED (local box e.g. 48 x 32 x 6144), and that is where blocking starts to matter.
  // Set -DSU2_NX/-DSU2_NY/-DSU2_NZ to reproduce a pencil on a single rank.
#if defined(SU2_NX) || defined(SU2_NY) || defined(SU2_NZ)
#ifndef SU2_NX
#define SU2_NX SU2_NGRID
#endif
#ifndef SU2_NY
#define SU2_NY SU2_NGRID
#endif
#ifndef SU2_NZ
#define SU2_NZ SU2_NGRID
#endif
  const device::IdxArray<NDim> boxSizes{{SU2_NX, SU2_NY, SU2_NZ}};
  const double sites = (double)SU2_NX * SU2_NY * SU2_NZ;
  auto toolBox = MemoryToolBox<NDim>::makeShared(boxSizes, nGhost, false);
#else
  const double sites = (double)nGrid * nGrid * nGrid;
  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost, false);
#endif
  toolBox->unsetVerbose();

  FieldCollection<SU2Field<T, NDim>, NDim, false, 1> U("U", toolBox);
  FieldCollection<SU2LieAlgebraField<T, NDim>, NDim, false, 1> Pi("Pi", toolBox);

  // Deterministic plane-wave IC, identical to benchmarks/su2_evolution.cpp (exact correctness gate).
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

  auto energy = [&]() {
    T elec = 2.0 * Total(i, 1, NDim, average(pow<2>(Pi(i)(1_c)) + pow<2>(Pi(i)(2_c)) + pow<2>(Pi(i)(3_c))));
    T mag = 2.0 * Total(i, 1, NDim, Total(j, 1, NDim, IfElse(i != j, average(1.0 - plaq(U, i, j).SU2Get(0_c)), 0.0)));
    return std::make_pair(elec, mag);
  };

  const auto [elec0, mag0] = energy();

#if KICK_VARIANT != 0
  // Hand-written site loop: the same assignment SU2LieAlgebraField::operator= performs, but with the
  // triple loop written out here instead of dispatched through Kokkos' MDRange policy. The
  // expression tree is built once per step and evaluated per site through the usual DoEval.
  auto assignManually = [&](auto &&target, auto &&r) {
    auto &f1 = target(1_c);
    auto &f2 = target(2_c);
    auto &f3 = target(3_c);

    f1.onBeforeAssignment(r); // config-space + ghost confirmation, exactly as operator= does
    f2.onBeforeAssignment(r);
    f3.onBeforeAssignment(r);
    PreGet::apply(r);

#if KICK_VARIANT == 2
#error "Variant 2 (raw-pointer leaves instead of Kokkos::View::operator()) was BUILT AND MEASURED, \
and it is a dead end: it matches the View to within noise. The View accessor is not the cost. The \
code was reverted rather than left behind a flag; see benchmarks/PERFORMANCE.md for the numbers and \
for what the remaining gap to the ceiling actually is (register spills -- try variants 4 and 5)."
#else
    const auto v1 = f1.getView();
    const auto v2 = f2.getView();
    const auto v3 = f3.getView();
#endif

    const auto layout = f1.getLayout();
    const auto Ls = layout.getLocalSizes();
    const device::Idx g = layout.getNGhosts();

    for (device::Idx i0 = g; i0 < g + (device::Idx)Ls[0]; ++i0)
      for (device::Idx i1 = g; i1 < g + (device::Idx)Ls[1]; ++i1)
#ifndef KICK_LADDER_NO_IVDEP
#pragma GCC ivdep
#endif
        for (device::Idx i2 = g; i2 < g + (device::Idx)Ls[2]; ++i2) {
          const auto result = DoEval::eval(r, i0, i1, i2);
          v1(i0, i1, i2) = result[1];
          v2(i0, i1, i2) = result[2];
          v3(i0, i1, i2) = result[3];
        }

    PostGet::apply(r);

    f1.setGhostsAreStale();
    f2.setGhostsAreStale();
    f3.setGhostsAreStale();
  };
#endif

  double kickSeconds = 0.0, driftSeconds = 0.0;
  for (size_t step = 0; step < nSteps; ++step) {
    auto t0 = std::chrono::steady_clock::now();
#if KICK_VARIANT == 4 || KICK_VARIANT == 5
    // Split the fused kick into smaller assignments. The fused form materializes all 12 four-link
    // products in one expression, whose live set does not fit in 32 zmm registers -- GCC spills, and
    // the kick executes 4.4x the loads of a hand-written kernel doing the same arithmetic. Splitting
    // trades a little extra Pi traffic (one read+write per assignment) for a smaller live set.
    //   4: one assignment per transverse direction j   (Pi -= dt (plaq - plaqBack))
    //   5: one assignment per term                     (Pi -= dt plaq;  Pi += dt plaqBack)
    for_in_range<1, NDim + 1>([&](auto i) {
      for_in_range<1, NDim + 1>([&](auto j) {
        if constexpr (decltype(i)::value != decltype(j)::value) {
#if KICK_VARIANT == 4
          Pi(i) = Pi(i) - dt * (plaq(U, i, j) - plaqBack(U, i, j));
#else
          Pi(i) = Pi(i) - dt * plaq(U, i, j);
          Pi(i) = Pi(i) + dt * plaqBack(U, i, j);
#endif
        }
      });
    });
#else
    for_in_range<1, NDim + 1>([&](auto i) {
      auto force = Pi(i) - dt * Total(j, 1, NDim, IfElse(i != j, plaq(U, i, j) - plaqBack(U, i, j), ZeroType()));
#if KICK_VARIANT == 0
      Pi(i) = force;
#else
      assignManually(Pi(i), force);
#endif
    });
#endif
    device::iteration::fence();
    auto t1 = std::chrono::steady_clock::now();
    for_in_range<1, NDim + 1>([&](auto i) { U(i) = exp(dt * Pi(i)) * U(i); });
    device::iteration::fence();
    auto t2 = std::chrono::steady_clock::now();

    kickSeconds += std::chrono::duration<double>(t1 - t0).count();
    driftSeconds += std::chrono::duration<double>(t2 - t1).count();
  }

  const auto [elec1, mag1] = energy();

  std::printf("variant  = %d\n", KICK_VARIANT);
  std::printf("N        = %zu\n", nGrid);
  std::printf("kick_ns  = %.1f\n", 1e9 * kickSeconds / nSteps / sites);
  std::printf("drift_ns = %.1f\n", 1e9 * driftSeconds / nSteps / sites);
  std::printf("plaq     = %.10g\n", (double)(mag1 / 12.0));
  std::printf("edrift   = %.3e\n", std::abs((elec1 + mag1) / (elec0 + mag0) - 1.0));
  return 0;
}

#endif
