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

  /**
   * @brief Exchange manager that routes ghost cell communication to P2P or MPI per (dimension, direction).
   *
   * On GPU builds with CUDA or HIP, the constructor probes which MPI neighbors reside on the same node
   * and have P2P-capable GPUs. For those neighbors, IPC handles are exchanged for the SEND buffers.
   * Each rank then READS from the remote's send buffer (pull model), matching ParaFaFT's proven approach.
   * Two MPI_Barrier calls on the shared-memory communicator synchronize the pack and read phases.
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
      mShmComm = shmComm;
      MPI_Comm_rank(mCartComm, &mMyRank);

      mP2PAvailable.fill(false);
      mRemoteSendUpPtr.fill(nullptr);
      mRemoteSendDownPtr.fill(nullptr);
      mRemoteHandleVersion.fill(0);

      probeP2P(shmComm);
#endif
    }

    ~ExchangeManager()
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      for (size_t i = 0; i < 2 * NDim; ++i) {
        if (mRemoteSendUpPtr[i] != nullptr) {
          p2p::ipcCloseHandle(mRemoteSendUpPtr[i]);
          mRemoteSendUpPtr[i] = nullptr;
        }
        if (mRemoteSendDownPtr[i] != nullptr) {
          p2p::ipcCloseHandle(mRemoteSendDownPtr[i]);
          mRemoteSendDownPtr[i] = nullptr;
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
    // Buffer handle exchange — call after (re)allocating send/recv buffers
    // ------------------------------------------------------------------

    void updateBufferHandles([[maybe_unused]] char *sendUpPtr, [[maybe_unused]] char *sendDownPtr,
                             [[maybe_unused]] char *recvUpPtr, [[maybe_unused]] char *recvDownPtr,
                             [[maybe_unused]] uint64_t version)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      exchangeIpcHandles(sendUpPtr, sendDownPtr, version);
#endif
    }

    // ------------------------------------------------------------------
    // Communication interface — ghost updater calls these
    // ------------------------------------------------------------------

    /**
     * @brief P2P pull-model exchange for one dimension.
     *
     * After packing + fence, call this instead of individual send/recv/waitall.
     * For P2P neighbors: signals "pack done" → remote reads from our send buffer →
     * remote signals "read done". For non-P2P neighbors: standard MPI Isend/Irecv.
     *
     * @param dimension     Which spatial dimension
     * @param sendUpPtr     Pointer to packed send-up slab (device memory)
     * @param sendDownPtr   Pointer to packed send-down slab (device memory)
     * @param recvUpPtr     Pointer to recv-up slab (device memory) — destination for pull
     * @param recvDownPtr   Pointer to recv-down slab (device memory) — destination for pull
     * @param byteCount     Number of bytes per slab
     * @param count         Number of MPI elements per slab
     * @param dataType      MPI datatype for non-P2P fallback
     */
    void exchange(size_t dimension, void *sendUpPtr, void *sendDownPtr, void *recvUpPtr, void *recvDownPtr,
                  size_t byteCount, int count, MPI_Datatype dataType)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      bool upP2P = isP2PUp(dimension);
      bool downP2P = isP2PDown(dimension);

      if (upP2P || downP2P) {
        // --- P2P pull model with barrier sync ---
        // Kokkos::fence was already called by the ghost updater (packing complete on GPU).

        // Post non-P2P MPI receives before the barrier (overlaps with barrier wait)
        if (!upP2P) mExchange.IrecvUp(dataType, dimension, recvUpPtr, count);
        if (!downP2P) mExchange.IrecvDown(dataType, dimension, recvDownPtr, count);

        // Barrier 1: all same-node ranks have finished packing → send buffers are safe to read
        MPI_Barrier(mShmComm);

        // Post non-P2P MPI sends (can overlap with P2P reads)
        if (!upP2P) mExchange.IsendUp(dataType, dimension, sendUpPtr, count);
        if (!downP2P) mExchange.IsendDown(dataType, dimension, sendDownPtr, count);

        // P2P reads from remote send buffers into local recv buffers
        if (upP2P)
          p2p::memcpyAsync(recvUpPtr, mRemoteSendUpPtr[dimension * 2 + 1], byteCount);
        if (downP2P)
          p2p::memcpyAsync(recvDownPtr, mRemoteSendDownPtr[dimension * 2 + 0], byteCount);

        // Ensure GPU reads complete
        p2p::streamSynchronize();

        // Wait for non-P2P MPI to complete
        if (!upP2P || !downP2P) mExchange.waitall();

        // Barrier 2: all same-node ranks have finished reading → send buffers safe to reuse
        MPI_Barrier(mShmComm);

        return;
      }
#endif
      // Pure MPI path (no P2P for this dimension)
      mExchange.IrecvUp(dataType, dimension, recvUpPtr, count);
      mExchange.IrecvDown(dataType, dimension, recvDownPtr, count);
      mExchange.IsendUp(dataType, dimension, sendUpPtr, count);
      mExchange.IsendDown(dataType, dimension, sendDownPtr, count);
      mExchange.waitall();
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
    MPI_Comm mShmComm = MPI_COMM_NULL;

    // Per (dimension, direction): indexed as [d * 2 + dir], dir: 0=up, 1=down
    std::array<bool, 2 * NDim> mP2PAvailable{};
    // IPC-mapped pointers to each neighbor's send buffers (we READ from these)
    // mRemoteSendUpPtr[d*2+dir]: the neighbor in direction 'dir' of dimension d's sendUp buffer
    std::array<void *, 2 * NDim> mRemoteSendUpPtr{};
    std::array<void *, 2 * NDim> mRemoteSendDownPtr{};
    std::array<uint64_t, 2 * NDim> mRemoteHandleVersion{};
    std::array<int, 2 * NDim> mNeighborRanks{};
    std::array<int, 2 * NDim> mNeighborDevices{};

    bool isP2PUp(size_t d) const { return mP2PAvailable[d * 2 + 0]; }
    bool isP2PDown(size_t d) const { return mP2PAvailable[d * 2 + 1]; }

    void probeP2P(MPI_Comm shmComm)
    {
      if (shmComm == MPI_COMM_NULL) return;

      MPI_Group worldGroup, shmGroup;
      MPI_Comm_group(mCartComm, &worldGroup);
      MPI_Comm_group(shmComm, &shmGroup);

      int shmSize;
      MPI_Comm_size(shmComm, &shmSize);

      std::vector<int> shmDevices(shmSize);
      MPI_Allgather(&mMyDevice, 1, MPI_INT, shmDevices.data(), 1, MPI_INT, shmComm);

      std::vector<int> shmGlobalRanks(shmSize);
      MPI_Allgather(&mMyRank, 1, MPI_INT, shmGlobalRanks.data(), 1, MPI_INT, shmComm);

      std::vector<std::pair<int, int>> rankDeviceMap;
      for (int i = 0; i < shmSize; ++i)
        rankDeviceMap.emplace_back(shmGlobalRanks[i], shmDevices[i]);

      auto &neighbours = mExchange.getNeighbours();

      for (size_t d = 0; d < NDim; ++d) {
        int upperNeighbor = neighbours.getUpperNeighbour(d);
        int lowerNeighbor = neighbours.getLowerNeighbour(d);

        mNeighborRanks[d * 2 + 0] = upperNeighbor;
        mNeighborRanks[d * 2 + 1] = lowerNeighbor;

        checkAndEnableP2P(d, 0, upperNeighbor, rankDeviceMap);
        checkAndEnableP2P(d, 1, lowerNeighbor, rankDeviceMap);
      }

      MPI_Group_free(&worldGroup);
      MPI_Group_free(&shmGroup);
    }

    void checkAndEnableP2P(size_t dim, int dir, int neighborRank,
                           const std::vector<std::pair<int, int>> &rankDeviceMap)
    {
      size_t idx = dim * 2 + dir;

      if (neighborRank == mMyRank) return;

      int neighborDevice = -1;
      for (auto &[rank, device] : rankDeviceMap) {
        if (rank == neighborRank) {
          neighborDevice = device;
          break;
        }
      }
      if (neighborDevice < 0) return;

      mNeighborDevices[idx] = neighborDevice;

      if (neighborDevice == mMyDevice) return;

      if (!p2p::canAccessPeer(mMyDevice, neighborDevice)) return;

      p2p::enablePeerAccess(neighborDevice);
      mP2PAvailable[idx] = true;

      sayMPI << "Ghost exchange: P2P enabled for dimension " << dim << (dir == 0 ? " (up)" : " (down)") << " to rank "
             << neighborRank << " (device " << neighborDevice << ")\n";
    }

    void exchangeIpcHandles(char *sendUpPtr, char *sendDownPtr, uint64_t version)
    {
      // Pull model: we need IPC handles for each neighbor's SEND buffers so we can READ from them.
      //
      // For "recvUp" (receiving data sent UP from our lower neighbor):
      //   - Our lower neighbor packed into their sendUp buffer
      //   - We need IPC handle for lower neighbor's sendUp buffer
      //   - We export our sendUp handle to our upper neighbor (they will recvUp = read our sendUp)
      //
      // For "recvDown" (receiving data sent DOWN from our upper neighbor):
      //   - Our upper neighbor packed into their sendDown buffer
      //   - We need IPC handle for upper neighbor's sendDown buffer
      //   - We export our sendDown handle to our lower neighbor (they will recvDown = read our sendDown)

      struct HandlePair {
        p2p::IpcHandlePacket sendUpPacket{};
        p2p::IpcHandlePacket sendDownPacket{};
        uint64_t version = 0;
      };

      for (size_t d = 0; d < NDim; ++d) {
        int upperRank = mNeighborRanks[d * 2 + 0];
        int lowerRank = mNeighborRanks[d * 2 + 1];

        // Exchange sendUp handles: we send ours to upper, receive lower's
        if (mP2PAvailable[d * 2 + 0] || mP2PAvailable[d * 2 + 1]) {
          // Pack our sendUp handle
          p2p::IpcHandlePacket mySendUpPacket{};
          if (sendUpPtr != nullptr) p2p::ipcGetHandle(sendUpPtr, mySendUpPacket.handle);
          mySendUpPacket.deviceId = mMyDevice;
          mySendUpPacket.version = version;

          // Pack our sendDown handle
          p2p::IpcHandlePacket mySendDownPacket{};
          if (sendDownPtr != nullptr) p2p::ipcGetHandle(sendDownPtr, mySendDownPacket.handle);
          mySendDownPacket.deviceId = mMyDevice;
          mySendDownPacket.version = version;

          // Exchange: send our sendUp handle to upper neighbor (they need it for their recvUp = read our sendUp)
          //           receive lower neighbor's sendUp handle (we need it for our recvUp = read their sendUp)
          p2p::IpcHandlePacket recvSendUpFromLower{};
          MPI_Status stat;
          int tag1 = 700 + d * 4 + 0;
          MPI_Sendrecv(&mySendUpPacket, sizeof(p2p::IpcHandlePacket), MPI_BYTE, upperRank, tag1, &recvSendUpFromLower,
                       sizeof(p2p::IpcHandlePacket), MPI_BYTE, lowerRank, tag1, mCartComm, &stat);

          // Exchange: send our sendDown handle to lower neighbor (they need it for their recvDown = read our sendDown)
          //           receive upper neighbor's sendDown handle (we need it for our recvDown = read their sendDown)
          p2p::IpcHandlePacket recvSendDownFromUpper{};
          int tag2 = 700 + d * 4 + 1;
          MPI_Sendrecv(&mySendDownPacket, sizeof(p2p::IpcHandlePacket), MPI_BYTE, lowerRank, tag2,
                       &recvSendDownFromUpper, sizeof(p2p::IpcHandlePacket), MPI_BYTE, upperRank, tag2, mCartComm,
                       &stat);

          // Open lower neighbor's sendUp handle (for our recvUp)
          if (mP2PAvailable[d * 2 + 0] && recvSendUpFromLower.version > mRemoteHandleVersion[d * 2 + 1]) {
            if (mRemoteSendUpPtr[d * 2 + 1] != nullptr) p2p::ipcCloseHandle(mRemoteSendUpPtr[d * 2 + 1]);
            mRemoteSendUpPtr[d * 2 + 1] =
                (recvSendUpFromLower.version > 0) ? p2p::ipcOpenHandle(recvSendUpFromLower.handle) : nullptr;
          }

          // Open upper neighbor's sendDown handle (for our recvDown)
          if (mP2PAvailable[d * 2 + 1] && recvSendDownFromUpper.version > mRemoteHandleVersion[d * 2 + 0]) {
            if (mRemoteSendDownPtr[d * 2 + 0] != nullptr) p2p::ipcCloseHandle(mRemoteSendDownPtr[d * 2 + 0]);
            mRemoteSendDownPtr[d * 2 + 0] =
                (recvSendDownFromUpper.version > 0) ? p2p::ipcOpenHandle(recvSendDownFromUpper.handle) : nullptr;
          }

          mRemoteHandleVersion[d * 2 + 0] = std::max(mRemoteHandleVersion[d * 2 + 0], recvSendDownFromUpper.version);
          mRemoteHandleVersion[d * 2 + 1] = std::max(mRemoteHandleVersion[d * 2 + 1], recvSendUpFromLower.version);
        }
      }
    }
#endif
  };

} // namespace TempLat::device_kokkos

#endif // HAVE_MPI

#endif
