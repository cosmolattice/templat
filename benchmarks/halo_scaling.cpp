/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

// ---------------------------------------------------------------------------
// Halo sub-phase diagnostic benchmark.
//
// Isolates the ghost/halo exchange from the leapfrog "kick" compute so the two
// terms can be compared independently across (N, ranks) on a controlled,
// FFT-free workload (pure real-space stencil; no Fourier transforms are ever
// executed, so this exercises only the halo path).
//
// Also measures the halo update of multi-component fields (ComplexField: 2
// components, SU2Field: 4 components) so the effect of coalescing the
// per-component exchanges into one message can be read off directly.
//
// Build:   cmake -B build-mpi -DTEMPLAT_BENCH=ON -DMPI=ON -DPARAFAFT=ON -DNOTHREADING=ON
//          cmake --build build-mpi --target bench-halo_scaling -j 8
// Run:     mpirun -n <ranks> ./build-mpi/benchmarks/bench-halo_scaling
//
// Grid overrides (compile-time, mirror the harness SIZES): -DHALO_NGRID=256
// -DHALO_NSTEPS=100 -DHALO_NGHOST=1 -DHALO_NDIM=3
// ---------------------------------------------------------------------------

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/session/sessionguard.h"
#include "TempLat/util/benchmark.h"
#include "TempLat/util/timer.h"

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfield.h"
#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/measuringtools/averager.h"

#ifndef HALO_NDIM
#define HALO_NDIM 3
#endif
#ifndef HALO_NGRID
#define HALO_NGRID 128
#endif
#ifndef HALO_NSTEPS
#define HALO_NSTEPS 100
#endif
#ifndef HALO_NGHOST
#define HALO_NGHOST 1
#endif

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = HALO_NDIM;
  using T = double;
  constexpr size_t nGrid = HALO_NGRID;
  constexpr size_t nGhost = HALO_NGHOST;
  constexpr size_t nSteps = HALO_NSTEPS;
  constexpr T dt = 0.01;

  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost, false);
  toolBox->unsetVerbose();

  int rank = 0, nRanks = 1;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
#endif

  Field<T, NDim> phi("phi", toolBox);
  Field<T, NDim> pi("pi", toolBox);
  ComplexField<T, NDim> psi("psi", toolBox);
  SU2Field<T, NDim> U("U", toolBox);

  // Seed with a non-trivial, communication-independent configuration.
  SpatialCoordinate x(toolBox);
  phi = x(1_c);
  pi = x(1_c);
  psi(0_c) = x(1_c);
  psi(1_c) = x(2_c);
  for_in_range<0, 4>([&](auto i) { U(i) = x(1_c); });

  // Bring phi's ghosts up to date once so the "kick" region below measures
  // pure stencil compute: phi is never written in the loop, so it stays fresh
  // and LatticeLaplacian(phi) triggers no exchange.
  phi.getMemoryManager()->confirmGhostsUpToDate();
  device::iteration::fence();

  // Report the normalization so ns/site/step can be derived from the CSV.
  auto layout = toolBox->mLayouts.getConfigSpaceLayout();
  size_t localSites = 1, globalSites = 1, interiorSites = 1;
  for (size_t d = 0; d < NDim; ++d) {
    const size_t ls = layout.getLocalSizes()[d];
    localSites *= ls;
    globalSites *= nGrid;
    interiorSites *= (ls > 2 * nGhost) ? (ls - 2 * nGhost) : 0; // sites whose stencil needs no ghost
  }
  const double interiorFraction = double(interiorSites) / double(localSites);
  auto decomp = toolBox->getDecomposition();
  if (rank == 0) {
    std::stringstream ss;
    ss << "halo_scaling: NDim=" << NDim << " nGrid=" << nGrid << " nGhost=" << nGhost << " nSteps=" << nSteps
       << " ranks=" << nRanks << "\n";
    ss << "  decomposition = [";
    for (size_t d = 0; d < NDim; ++d)
      ss << decomp[d] << (d + 1 < NDim ? "," : "");
    ss << "]  globalSites=" << globalSites << "  localSites/rank=" << localSites << "\n";
    ss << "  ns/site/step = Average[s] * 1e9 / localSites\n";
    sayMPI << ss.str();
  }

  Benchmark bench([&](Benchmark::Measurer &measurer) {
    for (size_t i = 0; i < nSteps; ++i) {
      // Pure stencil compute, ghosts already fresh -> no exchange inside.
      measurer.measure("kick", [&]() {
        pi = pi + dt * LatticeLaplacian(phi) * dt;
        device::iteration::fence();
      });
      // Isolated halo exchange: unconditional full ghost update of one field.
      measurer.measure("halo_update", [&]() {
        phi.updateGhosts();
        device::iteration::fence();
      });
      // Multi-component halo exchanges: 2 (complex) and 4 (SU2) components.
      // Divide by the component count to compare with halo_update.
      measurer.measure("cplx_halo", [&]() {
        psi.updateGhosts();
        device::iteration::fence();
      });
      measurer.measure("su2_halo", [&]() {
        U.updateGhosts();
        device::iteration::fence();
      });
    }
  });
  bench.run(1);
  bench.print();
  bench.log("halo_scaling");

  // -------------------------------------------------------------------------
  // Overlap ceiling (analytic). The interior of the kick (sites >= nGhost from
  // every local boundary) reads no ghost cells, so it could in principle run
  // WHILE the halo is in flight. kick_interior = kick * interiorFraction; if
  // the halo were perfectly overlapped behind that interior compute, the step
  // time would drop from (kick + halo) to (kick + halo - min(kick_interior,
  // halo)). This is only reported as a bound; no overlap is implemented.
  // -------------------------------------------------------------------------
  {
    const size_t nMeas = 30;
    // Warm + measure pure kick (phi ghosts stay fresh -> no exchange inside).
    for (size_t i = 0; i < 5; ++i) {
      pi = pi + dt * LatticeLaplacian(phi) * dt;
    }
    device::iteration::fence();
    Timer tk;
    for (size_t i = 0; i < nMeas; ++i)
      pi = pi + dt * LatticeLaplacian(phi) * dt;
    device::iteration::fence();
    const double kickNs = double(tk.nanoseconds()) / (nMeas * localSites);

    // Measure pure halo (unconditional full ghost update).
    for (size_t i = 0; i < 5; ++i)
      phi.updateGhosts();
    device::iteration::fence();
    Timer th;
    for (size_t i = 0; i < nMeas; ++i)
      phi.updateGhosts();
    device::iteration::fence();
    const double haloNs = double(th.nanoseconds()) / (nMeas * localSites);

    const double kickInterior = kickNs * interiorFraction;
    const double hidden = std::min(kickInterior, haloNs);
    const double stepNow = kickNs + haloNs;
    const double stepOverlap = stepNow - hidden;

    if (rank == 0) {
      std::stringstream ss;
      ss << "overlap ceiling, ns/site/step:\n";
      ss << "  kick              = " << kickNs << "\n";
      ss << "  halo              = " << haloNs << "\n";
      ss << "  interiorFraction  = " << interiorFraction << "\n";
      ss << "  kick_interior     = " << kickInterior << "\n";
      ss << "  step now          = " << stepNow << " (kick + halo)\n";
      ss << "  step w/ overlap   = " << stepOverlap << " (halo hidden behind interior compute)\n";
      ss << "  ceiling speedup   = " << stepNow / stepOverlap << "x\n";
      sayMPI << ss.str();
    }
  }
}
