#ifndef TEMPLAT_PARALLEL_KOKKOS_EXCHANGE_H
#define TEMPLAT_PARALLEL_KOKKOS_EXCHANGE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

#ifdef HAVE_MPI

#include "TempLat/parallel/mpi/cartesian/mpicartesianexchange.h"
#include "TempLat/parallel/mpi/mpitags.h"
#include "TempLat/util/log/saycomplete.h"

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#include "TempLat/parallel/devices/kokkos/kokkos_p2p.h"
#endif

#include <mpi.h>
#include <array>
#include <vector>

namespace TempLat::device_kokkos
{

  // MPI tag for P2P zero-byte completion signals — offset from ghost cell tag to avoid collisions
  namespace P2PTags
  {
    static constexpr int signalBase = 500;
  } // namespace P2PTags

  /**
   * @brief Exchange manager that routes ghost cell communication to P2P or MPI per (dimension, direction).
   *
   * On GPU builds with CUDA or HIP, the constructor probes which MPI neighbors reside on the same node
   * and have P2P-capable GPUs. For those neighbors, IPC handles are exchanged so that send operations
   * write directly into the remote rank's receive buffer. A zero-byte MPI message serves as the
   * completion signal. For all other neighbors, communication falls back to standard MPI_Isend/MPI_Irecv.
   *
   * On CPU builds, this class is a trivial wrapper around MPICartesianExchange with no overhead.
   */
  template <size_t NDim> class ExchangeManager
  {
  public:
    ExchangeManager(MPICartesianExchange exchange, [[maybe_unused]] MPI_Comm shmComm,
                    [[maybe_unused]] int myDeviceId)
        : mExchange(exchange)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      mMyDevice = myDeviceId;
      mCartComm = mExchange.getMPICartesianGroup().getComm();
      MPI_Comm_rank(mCartComm, &mMyRank);

      mP2PAvailable.fill(false);
      mRemoteRecvPtr.fill(nullptr);
      mRemoteHandleVersion.fill(0);

      probeP2P(shmComm);
#endif
    }

    ~ExchangeManager()
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      for (size_t i = 0; i < 2 * NDim; ++i) {
        if (mRemoteRecvPtr[i] != nullptr) {
          p2p::ipcCloseHandle(mRemoteRecvPtr[i]);
          mRemoteRecvPtr[i] = nullptr;
        }
      }
#endif
    }

    // Non-copyable (owns IPC handles)
    ExchangeManager(const ExchangeManager &) = delete;
    ExchangeManager &operator=(const ExchangeManager &) = delete;
    ExchangeManager(ExchangeManager &&) = default;
    ExchangeManager &operator=(ExchangeManager &&) = default;

    // ------------------------------------------------------------------
    // Buffer handle exchange — call after (re)allocating recv buffers
    // ------------------------------------------------------------------

    void updateBufferHandles([[maybe_unused]] char *recvUpPtr, [[maybe_unused]] char *recvDownPtr,
                             [[maybe_unused]] uint64_t version)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      exchangeIpcHandles(recvUpPtr, recvDownPtr, version);
#endif
    }

    // ------------------------------------------------------------------
    // Communication interface — ghost updater calls these
    // ------------------------------------------------------------------

    void recvUp(size_t dimension, void *ptr, int count, MPI_Datatype dataType)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      if (isP2PUp(dimension)) {
        // Post zero-byte signal receiver: the remote will signal after its P2P write completes
        char dummy = 0;
        MPI_Irecv(&dummy, 0, MPI_BYTE, mNeighborRanks[dimension * 2 + 0], signalTag(dimension, 0), mCartComm,
                   &mP2PSignalRecvRequests[dimension * 2 + 0]);
        return;
      }
#endif
      mExchange.IrecvUp(dataType, dimension, ptr, count);
    }

    void recvDown(size_t dimension, void *ptr, int count, MPI_Datatype dataType)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      if (isP2PDown(dimension)) {
        char dummy = 0;
        MPI_Irecv(&dummy, 0, MPI_BYTE, mNeighborRanks[dimension * 2 + 1], signalTag(dimension, 1), mCartComm,
                   &mP2PSignalRecvRequests[dimension * 2 + 1]);
        return;
      }
#endif
      mExchange.IrecvDown(dataType, dimension, ptr, count);
    }

    void sendUp(size_t dimension, void *ptr, size_t byteCount, [[maybe_unused]] int count,
                [[maybe_unused]] MPI_Datatype dataType)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      if (isP2PUp(dimension)) {
        // Direct write into remote rank's recvDown buffer (our up-neighbor's down-recv)
        p2p::memcpyAsync(mRemoteRecvPtr[dimension * 2 + 0], ptr, byteCount);
        return;
      }
#endif
      mExchange.IsendUp(dataType, dimension, ptr, count);
    }

    void sendDown(size_t dimension, void *ptr, size_t byteCount, [[maybe_unused]] int count,
                  [[maybe_unused]] MPI_Datatype dataType)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      if (isP2PDown(dimension)) {
        p2p::memcpyAsync(mRemoteRecvPtr[dimension * 2 + 1], ptr, byteCount);
        return;
      }
#endif
      mExchange.IsendDown(dataType, dimension, ptr, count);
    }

    void waitall([[maybe_unused]] size_t dimension)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      bool upP2P = isP2PUp(dimension);
      bool downP2P = isP2PDown(dimension);

      // Synchronize P2P: ensure async GPU copies are done, then signal remote
      if (upP2P || downP2P) p2p::streamSynchronize();

      if (upP2P) {
        char dummy = 0;
        MPI_Send(&dummy, 0, MPI_BYTE, mNeighborRanks[dimension * 2 + 0], signalTag(dimension, 0), mCartComm);
      }
      if (downP2P) {
        char dummy = 0;
        MPI_Send(&dummy, 0, MPI_BYTE, mNeighborRanks[dimension * 2 + 1], signalTag(dimension, 1), mCartComm);
      }

      // Wait for MPI transfers (non-P2P directions)
      if (!upP2P || !downP2P) {
        waitallMpiPartial(!upP2P, !downP2P);
      }

      // Wait for P2P signals from remote (confirms remote's write into our recv buffer is complete)
      if (upP2P) {
        MPI_Status stat;
        MPI_Wait(&mP2PSignalRecvRequests[dimension * 2 + 0], &stat);
      }
      if (downP2P) {
        MPI_Status stat;
        MPI_Wait(&mP2PSignalRecvRequests[dimension * 2 + 1], &stat);
      }
#else
      mExchange.waitall();
#endif
    }

    const MPICartesianGroup &getMPICartesianGroup() const { return mExchange.getMPICartesianGroup(); }

    // ------------------------------------------------------------------
    // Blocking exchange pass-through (used by CPU path)
    // ------------------------------------------------------------------

    void exchangeUp(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrSend, void *ptrReceive, int sendCount = 1)
    {
      mExchange.exchangeUp(dataType, dimension, ptrSend, ptrReceive, sendCount);
    }

    void exchangeDown(MPI_Datatype dataType, ptrdiff_t dimension, void *ptrSend, void *ptrReceive, int sendCount = 1)
    {
      mExchange.exchangeDown(dataType, dimension, ptrSend, ptrReceive, sendCount);
    }

  private:
    MPICartesianExchange mExchange;

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    int mMyDevice = -1;
    int mMyRank = -1;
    MPI_Comm mCartComm = MPI_COMM_NULL;

    // Per (dimension, direction): indexed as [d * 2 + dir], dir: 0=up, 1=down
    std::array<bool, 2 * NDim> mP2PAvailable{};
    std::array<void *, 2 * NDim> mRemoteRecvPtr{};
    std::array<uint64_t, 2 * NDim> mRemoteHandleVersion{};
    std::array<int, 2 * NDim> mNeighborRanks{};
    std::array<int, 2 * NDim> mNeighborDevices{};
    std::array<MPI_Request, 2 * NDim> mP2PSignalRecvRequests{};

    bool isP2PUp(size_t d) const { return mP2PAvailable[d * 2 + 0]; }
    bool isP2PDown(size_t d) const { return mP2PAvailable[d * 2 + 1]; }

    static int signalTag(size_t dimension, int direction) { return P2PTags::signalBase + dimension * 2 + direction; }

    void probeP2P(MPI_Comm shmComm)
    {
      if (shmComm == MPI_COMM_NULL) return;

      // Build a set of ranks in the shared-memory communicator
      MPI_Group worldGroup, shmGroup;
      MPI_Comm_group(mCartComm, &worldGroup);
      MPI_Comm_group(shmComm, &shmGroup);

      int shmSize;
      MPI_Comm_size(shmComm, &shmSize);

      // Allgather device IDs within the shared-memory communicator
      std::vector<int> shmDevices(shmSize);
      MPI_Allgather(&mMyDevice, 1, MPI_INT, shmDevices.data(), 1, MPI_INT, shmComm);

      // Allgather global ranks within shmComm so we can map shmRank -> globalRank -> device
      std::vector<int> shmGlobalRanks(shmSize);
      MPI_Allgather(&mMyRank, 1, MPI_INT, shmGlobalRanks.data(), 1, MPI_INT, shmComm);

      // Build globalRank -> deviceId map for same-node ranks
      std::vector<std::pair<int, int>> rankDeviceMap; // (globalRank, deviceId)
      for (int i = 0; i < shmSize; ++i)
        rankDeviceMap.emplace_back(shmGlobalRanks[i], shmDevices[i]);

      auto &neighbours = mExchange.getNeighbours();

      for (size_t d = 0; d < NDim; ++d) {
        int upperNeighbor = neighbours.getUpperNeighbour(d);
        int lowerNeighbor = neighbours.getLowerNeighbour(d);

        mNeighborRanks[d * 2 + 0] = upperNeighbor;
        mNeighborRanks[d * 2 + 1] = lowerNeighbor;

        // Check upper neighbor (up direction: we send to upper, recv from lower)
        // For P2P "up": we write into upper neighbor's recvDown buffer
        // The signal goes to/from the upper neighbor
        checkAndEnableP2P(d, 0, upperNeighbor, rankDeviceMap);

        // Check lower neighbor (down direction)
        checkAndEnableP2P(d, 1, lowerNeighbor, rankDeviceMap);
      }

      MPI_Group_free(&worldGroup);
      MPI_Group_free(&shmGroup);
    }

    void checkAndEnableP2P(size_t dim, int dir, int neighborRank,
                           const std::vector<std::pair<int, int>> &rankDeviceMap)
    {
      size_t idx = dim * 2 + dir;

      // Skip self (same rank means no MPI needed — handled by NOMPI path)
      if (neighborRank == mMyRank) return;

      // Check if neighbor is on the same node
      int neighborDevice = -1;
      for (auto &[rank, device] : rankDeviceMap) {
        if (rank == neighborRank) {
          neighborDevice = device;
          break;
        }
      }
      if (neighborDevice < 0) return; // not on same node

      mNeighborDevices[idx] = neighborDevice;

      // Check P2P capability
      if (neighborDevice == mMyDevice) return; // same device — no P2P needed (data is already accessible)

      if (!p2p::canAccessPeer(mMyDevice, neighborDevice)) return;

      p2p::enablePeerAccess(neighborDevice);
      mP2PAvailable[idx] = true;

      sayMPI << "Ghost exchange: P2P enabled for dimension " << dim << (dir == 0 ? " (up)" : " (down)")
             << " to rank " << neighborRank << " (device " << neighborDevice << ")\n";
    }

    void exchangeIpcHandles(char *recvUpPtr, char *recvDownPtr, uint64_t version)
    {
      // For each P2P-capable neighbor, exchange IPC handles for recv buffers.
      // Our upper neighbor needs a handle to our recvUp buffer (they write into it as their sendDown).
      // Wait — the directions are:
      //   - We sendUp to upper neighbor. Upper neighbor's recvDown is where our data lands.
      //     So we need IPC handle for upper neighbor's recvDown buffer.
      //   - Upper neighbor sendDown to us. Our recvUp is where their data lands.
      //     So upper neighbor needs IPC handle for our recvUp buffer.
      //
      // In our "push" model:
      //   - sendUp writes to remote's recv buffer. We need the remote's recv handle.
      //   - The remote's recv for our sendUp is: the remote's recvDown (since we send up, they receive down).
      //     But from the remote's perspective, when we sendUp, the remote receives from its lower neighbor.
      //     The remote posts recvUp from its lower neighbor... wait, let me think about this clearly.
      //
      // Ghost exchange semantics in MPICartesianExchange:
      //   - IsendUp sends to getUpperNeighbour, IrecvUp recvs from getLowerNeighbour
      //   - IsendDown sends to getLowerNeighbour, IrecvDown recvs from getUpperNeighbour
      //
      // So when rank A sendUp to rank B (B = A's upper neighbor):
      //   - A is B's lower neighbor
      //   - B will IrecvUp from its lower neighbor = A, into B's recvUpSlab
      //
      // For P2P: A needs to write directly into B's recvUpSlab.
      // So A needs the IPC handle for B's recvUp buffer.
      //
      // Similarly, when A sendDown to rank C (C = A's lower neighbor):
      //   - A is C's upper neighbor
      //   - C will IrecvDown from its upper neighbor = A, into C's recvDownSlab
      //   - So A needs the IPC handle for C's recvDown buffer.
      //
      // Exchange protocol:
      //   - Each rank sends its recvUp handle to its lower neighbor (who will sendUp to us)
      //   - Each rank sends its recvDown handle to its upper neighbor (who will sendDown to us)

      for (size_t d = 0; d < NDim; ++d) {
        // --- Up direction: we sendUp to upper neighbor, need their recvUp handle ---
        if (mP2PAvailable[d * 2 + 0]) {
          int upperRank = mNeighborRanks[d * 2 + 0];

          // Send our recvUp handle to our lower neighbor (they sendUp to us, need our recvUp)
          // Receive upper neighbor's recvUp handle (we sendUp to them, need their recvUp)
          p2p::IpcHandlePacket sendPacket{};
          if (recvUpPtr != nullptr) p2p::ipcGetHandle(recvUpPtr, sendPacket.handle);
          sendPacket.deviceId = mMyDevice;
          sendPacket.version = version;

          p2p::IpcHandlePacket recvPacket{};
          MPI_Status stat;
          // We send our recvUp to lower neighbor, receive upper neighbor's recvUp
          // But actually: upper neighbor sends us their recvUp handle
          // Use MPI_Sendrecv: send our recvUp handle to lower, recv upper's recvUp handle from upper
          int lowerRank = mNeighborRanks[d * 2 + 1]; // our lower neighbor
          // Lower neighbor will sendUp to us and needs our recvUp handle
          // Upper neighbor: we sendUp to them and need their recvUp handle (they send it to us)
          MPI_Sendrecv(&sendPacket, sizeof(p2p::IpcHandlePacket), MPI_BYTE, lowerRank, signalTag(d, 0) + 100,
                       &recvPacket, sizeof(p2p::IpcHandlePacket), MPI_BYTE, upperRank, signalTag(d, 0) + 100,
                       mCartComm, &stat);

          if (recvPacket.version > mRemoteHandleVersion[d * 2 + 0]) {
            if (mRemoteRecvPtr[d * 2 + 0] != nullptr) p2p::ipcCloseHandle(mRemoteRecvPtr[d * 2 + 0]);
            mRemoteRecvPtr[d * 2 + 0] = (recvPacket.version > 0) ? p2p::ipcOpenHandle(recvPacket.handle) : nullptr;
            mRemoteHandleVersion[d * 2 + 0] = recvPacket.version;
          }
        }

        // --- Down direction: we sendDown to lower neighbor, need their recvDown handle ---
        if (mP2PAvailable[d * 2 + 1]) {
          int lowerRank = mNeighborRanks[d * 2 + 1];

          p2p::IpcHandlePacket sendPacket{};
          if (recvDownPtr != nullptr) p2p::ipcGetHandle(recvDownPtr, sendPacket.handle);
          sendPacket.deviceId = mMyDevice;
          sendPacket.version = version;

          p2p::IpcHandlePacket recvPacket{};
          MPI_Status stat;
          int upperRank = mNeighborRanks[d * 2 + 0]; // our upper neighbor
          // Upper neighbor will sendDown to us and needs our recvDown handle
          // Lower neighbor: we sendDown to them and need their recvDown handle
          MPI_Sendrecv(&sendPacket, sizeof(p2p::IpcHandlePacket), MPI_BYTE, upperRank, signalTag(d, 1) + 100,
                       &recvPacket, sizeof(p2p::IpcHandlePacket), MPI_BYTE, lowerRank, signalTag(d, 1) + 100,
                       mCartComm, &stat);

          if (recvPacket.version > mRemoteHandleVersion[d * 2 + 1]) {
            if (mRemoteRecvPtr[d * 2 + 1] != nullptr) p2p::ipcCloseHandle(mRemoteRecvPtr[d * 2 + 1]);
            mRemoteRecvPtr[d * 2 + 1] = (recvPacket.version > 0) ? p2p::ipcOpenHandle(recvPacket.handle) : nullptr;
            mRemoteHandleVersion[d * 2 + 1] = recvPacket.version;
          }
        }
      }
    }

    void waitallMpiPartial(bool waitUp, bool waitDown)
    {
      // Collect active MPI requests and wait on them
      // mExchange uses indices: 0=IrecvUp, 1=IrecvDown, 2=IsendUp, 3=IsendDown
      // We just call the full waitall — inactive requests are MPI_REQUEST_NULL and are skipped
      // However, we can only call waitall if we actually posted MPI requests.
      // For safety, just call the full waitall which handles MPI_REQUEST_NULL correctly.
      if (waitUp || waitDown) mExchange.waitall();
    }
#endif
  };

} // namespace TempLat::device_kokkos

#endif // HAVE_MPI

#endif
