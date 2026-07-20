
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesiangroup.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"
#include "TempLat/util/tdd/tdd.h"

#include <vector>

namespace TempLat
{
  /** @brief Tests that the Cartesian group adopts the FFT backend's topology rather than
   *  inventing its own.
   *
   *  The bug these guard against: TempLat and the FFT backend each called MPI_Cart_create with
   *  reorder = 1 and reconciled only their grid *shapes*. Local starts came from the backend's
   *  coordinates, ghost neighbours from MPI_Cart_shift on TempLat's communicator, and nothing
   *  forced the two rank->coordinate mappings to agree. It worked only because mainstream MPI
   *  implementations ignore `reorder`; when they disagree the failure is silent wrong data at
   *  subdomain boundaries.
   */
  struct MPICartesianGroupTester {

    /** getPosition() must be exactly what MPI_Cart_coords reports on getComm(). */
    template <size_t NDim> static bool positionMatchesCartCoords(const device::IdxArray<NDim> &nGrid)
    {
      auto group = FFTMPIDomainSplit<NDim>::makeMPIGroup(MPICommReference(), nGrid);

      std::vector<int> coords(NDim, -1);
#ifdef HAVE_MPI
      int rank = 0;
      MPI_Comm_rank(group.getComm(), &rank);
      if (MPI_Cart_coords(group.getComm(), rank, static_cast<int>(NDim), coords.data()) != 0) return false;
#else
      std::fill(coords.begin(), coords.end(), 0);
#endif

      const auto &position = group.getPosition();
      if (position.size() != NDim) return false;
      for (size_t i = 0; i < NDim; ++i)
        if (position[i] != coords[i]) return false;
      return true;
    }

    /** With reorder = false the Cartesian communicator must preserve the base rank order.
     *
     *  If this ever fails, every consumer that mixes the two rank spaces is wrong: getRank()
     *  returns the Cartesian rank while MemoryToolBox::getMPIRank()/amIRoot() return the base one.
     */
    template <size_t NDim> static bool cartesianRankEqualsBaseRank(const device::IdxArray<NDim> &nGrid)
    {
      auto group = FFTMPIDomainSplit<NDim>::makeMPIGroup(MPICommReference(), nGrid);
#ifdef HAVE_MPI
      int cartRank = 0;
      MPI_Comm_rank(group.getComm(), &cartRank);
      return cartRank == group.getBaseComm().rank();
#else
      (void)group;
      return true;
#endif
    }

    /** The group's coordinates must equal the coordinates the FFT backend reports for this rank.
     *
     *  This is the actual invariant: the layout's local starts are derived from the backend's
     *  coordinates, so if these differ the group services the wrong subdomain.
     */
    template <size_t NDim> static bool coordsMatchBackendTopology(const device::IdxArray<NDim> &nGrid)
    {
      MPICommReference baseComm;
      const auto topology = FFTLibrarySelector<NDim>::topology(baseComm, nGrid);
      auto group = FFTMPIDomainSplit<NDim>::makeMPIGroup(baseComm, nGrid);

      const auto &position = group.getPosition();
      for (size_t i = 0; i < NDim; ++i) {
        if (position[i] != topology.coords[i]) return false;
        if (group.getDecomposition()[i] != topology.dims[i]) return false;
      }
      return true;
    }

    /** The coordinate cross-check must actually fire when the topology disagrees.
     *
     *  No mainstream MPI honours `reorder`, so the original bug cannot reproduce naturally and
     *  the guard would otherwise never be exercised. Build a communicator whose rank order is
     *  REVERSED relative to the world comm, then hand the group the coordinates computed for the
     *  un-reversed order. That is precisely the "two differently-ordered communicators over the
     *  same ranks" situation, and it must throw rather than silently proceed.
     */
    template <size_t NDim> static bool rejectsDisagreeingCoords(const device::IdxArray<NDim> &nGrid)
    {
#ifdef HAVE_MPI
      MPICommReference world;
      const int size = static_cast<int>(world.size());
      if (size < 2) return true; // reversal is a no-op at one rank

      const auto topology = FFTLibrarySelector<NDim>::topology(world, nGrid);

      // Reverse the rank order: rank r in world becomes rank (size-1-r) here.
      MPI_Comm reversedRaw;
      MPI_Comm_split(world, 0, size - 1 - world.rank(), &reversedRaw);
      // Hand ownership to MPICommReference and never free it by hand: it is refcounted and
      // releases the communicator itself. Doing both is a double free.
      MPICommReference reversed(reversedRaw);

      // Coordinates as computed on the ORIGINAL ordering, applied to the reversed communicator.
      std::vector<int> staleCoords(topology.coords.begin(), topology.coords.end());

      bool threw = false;
      try {
        MPICartesianGroup group(reversed, static_cast<device::Idx>(NDim), topology.dimsVector(), staleCoords);
      } catch (const MPICartesianGroupException &) {
        threw = true;
      } catch (...) {
        return false; // wrong exception type
      }

      // A rank whose reversed position happens to coincide with its original one legitimately
      // does not throw (the middle rank of an odd count, for instance). Only require that the
      // guard fired somewhere in the job — otherwise it is not actually wired up.
      int localThrew = threw ? 1 : 0;
      int anyThrew = 0;
      MPI_Allreduce(&localThrew, &anyThrew, 1, MPI_INT, MPI_SUM, world);
      return anyThrew > 0;
#else
      (void)nGrid;
      return true;
#endif
    }

    /** A decomposition that does not tile the communicator must be rejected outright.
     *
     *  This used to be "repaired" by overwriting entry 0 and zeroing trailing entries from the
     *  right, which silently turned an explicitly-pinned backend grid into a different one.
     */
    static bool rejectsNonTilingDecomposition()
    {
#ifdef HAVE_MPI
      MPICommReference world;
      if (world.size() < 2) return true;

      std::vector<int> tooSmall(3, 1); // product 1, but the communicator has >1 rank
      try {
        MPICartesianGroup group(world, 3, tooSmall);
      } catch (const MPICartesianGroupException &) {
        return true;
      } catch (...) {
        return false;
      }
      return false; // must not silently "repair"
#else
      return true;
#endif
    }

    template <size_t NDim, class... Vs> static void runCase(TDDAssertion &tdd, Vs... vs)
    {
      device::IdxArray<NDim> nGrid{static_cast<device::Idx>(vs)...};
      tdd.verify(positionMatchesCartCoords<NDim>(nGrid));
      tdd.verify(cartesianRankEqualsBaseRank<NDim>(nGrid));
      tdd.verify(coordsMatchBackendTopology<NDim>(nGrid));
      tdd.verify(rejectsDisagreeingCoords<NDim>(nGrid));
    }

    static void Test(TDDAssertion &tdd)
    {
      runCase<2>(tdd, 16, 16);
      runCase<3>(tdd, 16, 16, 16);
      runCase<3>(tdd, 8, 32, 16); // anisotropic: exercises the slab-vs-pencil heuristic
      tdd.verify(rejectsNonTilingDecomposition());
    }
  };

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MPICartesianGroupTester> test;
}
