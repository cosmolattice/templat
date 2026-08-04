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

  // MemoryToolBox manages the memory layout and ghost cells for the fields.
  auto toolBox = MemoryToolBox<3>::makeShared(128, 1);
  toolBox->unsetVerbose();

  // Create a random Gaussian field phi with mean 2 and standard deviation 1.
  Field<double, 3> phi("phi", toolBox);
  phi = 2 + RandomGaussianFieldConfig<double, 3>("seed", toolBox);

  // Compute the average of phi across the lattice.
  const double result = average(phi);

  int rank = 0;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == 0) sayMPI << "Hello, TempLat!  <phi> = " << result << "\n";

  return 0;
};
