/* phi^4 scalar field leapfrog evolution.
 *
 * Per step (kick-drift leapfrog, dx = 1):
 *   kick : pi  = pi  + dt * ( lap(phi) - m2*phi - lambda*phi^3 )
 *   drift: phi = phi + dt * pi
 *
 * Usage: classical_phi4 [N] [steps] (defaults 64 100)
 *
 * Output: CSV lines
 *   t,phi_avg,phi2_avg,pi_avg,pi2_avg
 */

#include "TempLat/session/sessionguard.h"

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/measuringtools/averager.h"

#include <cstdio>
#include <cstdlib>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
  using namespace TempLat;

  // The SessionGuard initializes MPI and Kokkos, and finalizes them at the end of main().
  SessionGuard guard(argc, argv, false);

  // Compile-time constants: float type and number of dimensions.
  using T = double;
  constexpr size_t NDim = 3;

  // Runtime parameters: grid size and number of steps.
  const size_t nGrid = argc > 1 ? std::atol(argv[1]) : 64;
  const int nSteps = argc > 2 ? std::atoi(argv[2]) : 100;
  constexpr size_t nGhost = 1;

  const T dt = 0.02;
  const T m2 = 1.0;
  const T lambda = 0.1;

  // MemoryToolBox manages the memory layout and ghost cells for the fields.
  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost, false);
  toolBox->unsetVerbose();

  // Fields for the scalar field phi and its conjugate momentum pi.
  Field<T, NDim> phi("phi", toolBox);
  Field<T, NDim> pi("pi", toolBox);

  // Random initial condition
  phi.inFourierSpace() = RandomGaussianField<T, NDim>("bench_phi4_seed", toolBox);
  // Add a zero mode to phi to avoid starting from a pure vacuum configuration
  // Normalization is from Fourier transform conventions: phi(x) = sum_k phi(k) exp(ikx), so the zero mode in Fourier
  // space is multiplied by the volume when transforming back to configuration space.
  phi.inFourierSpace().setZeroMode(pow<NDim>(nGrid) * 1.0);
  // Initialize pi to zero
  pi = T(0);

  // The two phases of the leapfrog evolution: kick, drift
  auto kick = [&]() { pi = pi + dt * (LatticeLaplacian(phi) - m2 * phi - lambda * pow<3>(phi)); };
  auto drift = [&]() { phi = phi + dt * pi; };

  std::vector<double> av_phi, av_phi2, av_pi, av_pi2;

  // Perform the time evolution for nSteps
  for (int i = 0; i < nSteps; ++i) {
    // Compute averages every 10 steps
    if (i % 10 == 0) {
      av_phi.push_back(average(phi));
      av_phi2.push_back(average(pow<2>(phi)));
      av_pi.push_back(average(pi));
      av_pi2.push_back(average(pow<2>(pi)));
    }

    kick();
    drift();
  }

  int rank = 0;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  // on rank 0, print the time series of averages
  if (rank == 0) {
    std::printf("t,phi_avg,phi2_avg,pi_avg,pi2_avg\n");
    for (size_t i = 0; i < av_phi.size(); ++i)
      std::printf("%zu,%g,%g,%g,%g\n", i * 10, av_phi[i], av_phi2[i], av_pi[i], av_pi2[i]);
  }

  return 0;
}
