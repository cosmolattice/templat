#ifndef TEMPLAT_PARALLEL_MPI_CARTESIAN_MPICARTESIANEXCHANGE_H
#define TEMPLAT_PARALLEL_MPI_CARTESIAN_MPICARTESIANEXCHANGE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019

#include "TempLat/parallel/mpi/cartesian/mpicartesiangroup.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesianneighbours.h"
#include "TempLat/parallel/mpi/mpitags.h"

namespace TempLat
{

  /** @brief A class which handles the exchange between neighbours in the cartesian group.
   *  Has two methods: exchangeUp and exchangeDown, which take as input a datatype, a dimension
   *  (because you need to specify what is up and down),
   *  a pointer for the sending memory and a pointer
   *  for the receiving memory.
   *
   * Unit test: ctest -R test-mpicartesianexchange
   **/
  class MPICartesianExchange
  {
  public:
    // Put public methods here. These should change very little over time.
    MPICartesianExchange(MPICartesianGroup group) : mGroup(group), mNeighbours(mGroup)
    {
#ifdef HAVE_MPI
      // The P2P-mixed exchange path may post only a subset of {IsendUp, IrecvUp, IsendDown, IrecvDown}.
      // waitall() always passes all 4 slots to MPI_Waitall, which only tolerates unused slots if they
      // are MPI_REQUEST_NULL (not uninitialized memory).
      mRequests.fill(MPI_REQUEST_NULL);
#endif
    }

    void exchangeUp(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrSend, void *ptrReceive, int sendCount = 1)
    {
#ifdef HAVE_MPI
      MPI_Status stat;
      MPI_Sendrecv(ptrSend, sendCount, dataType, mNeighbours.getUpperNeighbour(dimension), MPITags::dataShiftGhostCells,
                   ptrReceive, sendCount, dataType, mNeighbours.getLowerNeighbour(dimension),
                   MPITags::dataShiftGhostCells, mGroup.getComm(), &stat);
#endif
    }

    void exchangeDown(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrSend, void *ptrReceive, int sendCount = 1)
    {
#ifdef HAVE_MPI
      MPI_Status stat;
      MPI_Sendrecv(ptrSend, sendCount, dataType, mNeighbours.getLowerNeighbour(dimension), MPITags::dataShiftGhostCells,
                   ptrReceive, sendCount, dataType, mNeighbours.getUpperNeighbour(dimension),
                   MPITags::dataShiftGhostCells, mGroup.getComm(), &stat);
#endif
    }

    void IsendUp(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrSend, int sendCount = 1)
    {
#ifdef HAVE_MPI
      MPI_Isend(ptrSend, sendCount, dataType, mNeighbours.getUpperNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[2]);
#endif
    }
    void IrecvUp(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrReceive, int sendCount = 1)
    {
#ifdef HAVE_MPI
      MPI_Irecv(ptrReceive, sendCount, dataType, mNeighbours.getLowerNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[0]);
#endif
    }
    void IsendDown(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrSend, int sendCount = 1)
    {
#ifdef HAVE_MPI
      MPI_Isend(ptrSend, sendCount, dataType, mNeighbours.getLowerNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[3]);
#endif
    }
    void IrecvDown(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrReceive, int sendCount = 1)
    {
#ifdef HAVE_MPI
      MPI_Irecv(ptrReceive, sendCount, dataType, mNeighbours.getUpperNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[1]);
#endif
    }
    void waitall()
    {
#ifdef HAVE_MPI
      std::array<MPI_Status, 4> stats;
      MPI_Waitall(4, mRequests.data(), stats.data());
#endif
    }

#ifdef HAVE_MPI
    /** @brief Non-blocking up/down exchange where each of the four messages carries its own absolute-address
     *  datatype (built by the caller, typically MPI_Type_create_hindexed_block over several component faces)
     *  and is addressed with MPI_BOTTOM. Posting order (recvs before sends, Up before Down) matches
     *  waitall()/exchangeUpDownNonBlocking so the single ghost tag stays matched on two-rank dimensions. */
    void exchangeUpDownBottom(ptrdiff_t dimension, MPI_Datatype sendUpType, MPI_Datatype recvUpType,
                              MPI_Datatype sendDownType, MPI_Datatype recvDownType)
    {
      MPI_Irecv(MPI_BOTTOM, 1, recvUpType, mNeighbours.getLowerNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[0]);
      MPI_Irecv(MPI_BOTTOM, 1, recvDownType, mNeighbours.getUpperNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[1]);
      MPI_Isend(MPI_BOTTOM, 1, sendUpType, mNeighbours.getUpperNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[2]);
      MPI_Isend(MPI_BOTTOM, 1, sendDownType, mNeighbours.getLowerNeighbour(dimension), MPITags::dataShiftGhostCells,
                mGroup.getComm(), &mRequests[3]);
      waitall();
    }
#endif

    const MPICartesianGroup &getMPICartesianGroup() const { return mGroup; }

    MPICartesianNeighbours &getNeighbours() { return mNeighbours; }
    const MPICartesianNeighbours &getNeighbours() const { return mNeighbours; }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    MPICartesianGroup mGroup;
    MPICartesianNeighbours mNeighbours;
#ifdef HAVE_MPI
    std::array<MPI_Request, 4> mRequests;
#endif
  };
} // namespace TempLat

#endif
