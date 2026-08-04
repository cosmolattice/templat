/* SU(2) classical Yang-Mills leapfrog evolution.
 *
 * Pure-gauge Hamiltonian evolution in temporal gauge (A_0 = 0), dx = 1,
 * kick-drift leapfrog. Per step:
 *   kick : pi_i -= dt * sum_{j != i} (plaq(U,i,j) - plaqBack(U,i,j))
 *          (assignment to the Lie-algebra field drops the identity component,
 *           i.e. projects onto the algebra)
 *   drift: U_i = exp(dt * pi_i) * U_i          (exact on the group)
 *
 *  The algebra fields store p with U-dot = (i p.sigma) U;
 *
 * Initial condition (deterministic):
 *   U_1 = exp(i A sin(2 pi x_1 / N) sigma_1)   (angle varies along dir 2)
 *   U_2 = exp(i A sin(2 pi x_2 / N) sigma_2)   (angle varies along dir 3)
 *   U_3 = exp(i A sin(2 pi x_0 / N) sigma_3)   (angle varies along dir 1)
 * with A = 0.5, pi = 0. The plane waves have period exactly N.
 *
 * Conserved energy per site (E = i p.sigma antihermitean):
 *   H = 2 sum_{i,a} p_ia^2 + 2 sum_{i != j} (1 - c0(P_ij))
 * where c0 is the identity component, c0(P) = Re Tr P / 2.
 *
 * Usage: classical_su2 [N] [steps]   (defaults 64 100)
 *
 * Output: CSV lines
 *   t,elec_avg,mag_avg,edrift
 */

#include "TempLat/session/sessionguard.h"

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

  const T dt = 0.001;
  const T amp = 0.5;

  // MemoryToolBox manages the memory layout and ghost cells for the fields.
  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost, false);
  toolBox->unsetVerbose();

  // Links (4 real components per direction) and algebra momenta (3 per direction).
  FieldCollection<SU2Field<T, NDim>, NDim, false, 1> U("U", toolBox);
  FieldCollection<SU2LieAlgebraField<T, NDim>, NDim, false, 1> Pi("Pi", toolBox);

  SpatialCoordinate<NDim> x(toolBox);
  const T w = 2.0 * M_PI / nGrid;

  // U_1 = exp(i th sigma_1), th = amp*sin(w*x_1): (c0,c1,c2,c3) = (cos th, sin th, 0, 0)
  U(1_c)(0_c) = cos(amp * sin(w * getVectorComponent(x, Tag<1>())));
  U(1_c)(1_c) = sin(amp * sin(w * getVectorComponent(x, Tag<1>())));
  U(1_c)(2_c) = T(0);
  U(1_c)(3_c) = T(0);
  // U_2 = exp(i th sigma_2), th = amp*sin(w*x_2)
  U(2_c)(0_c) = cos(amp * sin(w * getVectorComponent(x, Tag<2>())));
  U(2_c)(1_c) = T(0);
  U(2_c)(2_c) = sin(amp * sin(w * getVectorComponent(x, Tag<2>())));
  U(2_c)(3_c) = T(0);
  // U_3 = exp(i th sigma_3), th = amp*sin(w*x_0)
  U(3_c)(0_c) = cos(amp * sin(w * getVectorComponent(x, Tag<0>())));
  U(3_c)(1_c) = T(0);
  U(3_c)(2_c) = T(0);
  U(3_c)(3_c) = sin(amp * sin(w * getVectorComponent(x, Tag<0>())));

  // Pi = 0
  for_in_range<1, NDim + 1>([&](auto i) {
    Pi(i)(1_c) = T(0);
    Pi(i)(2_c) = T(0);
    Pi(i)(3_c) = T(0);
  });

  // Per-site energy H = elec + mag (see header); conserved up to O(dt^2).
  auto energy_per_site = [&]() {
    T elec = 2.0 * Total(i, 1, NDim, average(pow<2>(Pi(i)(1_c)) + pow<2>(Pi(i)(2_c)) + pow<2>(Pi(i)(3_c))));
    T mag = 2.0 * Total(i, 1, NDim, Total(j, 1, NDim, IfElse(i != j, average(1.0 - plaq(U, i, j).SU2Get(0_c)), 0.0)));
    return std::make_pair(elec, mag);
  };
  const auto [elec0, mag0] = energy_per_site();
  const T energy0 = elec0 + mag0;

  // The two phases of the leapfrog evolution: kick, drift
  auto kick = [&]() {
    for_in_range<1, NDim + 1>([&](auto i) {
      Pi(i) = Pi(i) - dt * Total(j, 1, NDim, IfElse(i != j, plaq(U, i, j) - plaqBack(U, i, j), ZeroType()));
    });
  };
  auto drift = [&]() { for_in_range<1, NDim + 1>([&](auto i) { U(i) = exp(dt * Pi(i)) * U(i); }); };

  std::vector<double> av_elec, av_mag, av_edrift;

  for (int i = 0; i < nSteps; ++i) {
    // Compute averages every 10 steps
    if (i % 10 == 0) {
      const auto [elec1, mag1] = energy_per_site();
      const T edrift = std::abs((elec1 + mag1) / energy0 - 1.0);

      av_elec.push_back(elec1);
      av_mag.push_back(mag1);
      av_edrift.push_back(edrift);
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
    std::printf("t,elec_avg,mag_avg,edrift\n");
    for (size_t i = 0; i < av_elec.size(); ++i)
      std::printf("%zu,%g,%g,%g\n", i * 10, av_elec[i], av_mag[i], av_edrift[i]);
  }

  return 0;
}
