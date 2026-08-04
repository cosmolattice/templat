/* Langevin dynamics example: solves the stochastic PDE
 * d phi / dt = Laplacian(phi) - lambda/2 phi^3 + eta
 * where eta is a Gaussian noise with <eta(x,t) eta(x',t')> = delta(x-x') delta(t-t').
 *
 * Usage: langevin [N] [steps]   (defaults 128 100)
 *
 * Output: average of phi every 10 steps
 */

#include "TempLat.h"

int main(int argc, char *argv[])
{
  using namespace TempLat;

  // The SessionGuard initializes MPI and Kokkos, and finalizes them at the end of main().
  SessionGuard guard(argc, argv);

  // Compile-time constants: float type and number of dimensions.
  using T = double;
  constexpr size_t NDim = 3;

  // Runtime parameters: grid size.
  const size_t nGrid = argc > 1 ? std::atol(argv[1]) : 128;
  const int nSteps = argc > 2 ? std::atoi(argv[2]) : 100;
  const size_t ghosts = 1;

  const T lambda = 0.1;
  const T dt = 0.01;

  // MemoryToolBox manages the memory layout and ghost cells for the fields.
  auto toolBox = MemoryToolBox<3>::makeShared(nGrid, ghosts, false);
  toolBox->unsetVerbose();

  // Field classes hold data
  Field<T, NDim> phi("phi", toolBox);

  // Create an object giving gaussian fluctuations
  RandomGaussianFieldConfig<T, NDim> eta("seed", toolBox);

  // Create an abstract expression that can be  evaluated on the lattice; the noise is rescaled
  // by 1/sqrt(dt) to obtain the discretized delta correlation in time
  auto expr = LatLapl(phi) - lambda / 2 * pow<3>(phi) + device::sqrt(1 / dt) * eta;

  int rank = 0;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  // Now perform 100 evolution steps
  for (int step = 0; step < nSteps; ++step) {
    // every 10 steps, print the average of phi
    if (step % 10 == 0) {
      auto result = average(phi);
      if (rank == 0) std::printf("Step %zu, <phi> = %g\n", step, result);
    }
    // This evaluates the expression and assigns it to the data held by phi
    phi = phi + dt * expr;
  }

  return 0;
};
