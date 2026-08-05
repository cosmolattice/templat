#ifndef TEMPLAT_PARALLEL_KOKKOS_EXCHANGE_H
#define TEMPLAT_PARALLEL_KOKKOS_EXCHANGE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#ifdef HAVE_MPI

#include "TempLat/parallel/mpi/cartesian/mpicartesianexchange.h"
#include "TempLat/util/log/saycomplete.h"

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#include "TempLat/parallel/devices/kokkos/kokkos_p2p.h"
#endif

#include <mpi.h>
#include <algorithm>
#include <array>
#include <cstring>
#include <exception>
#include <vector>

namespace TempLat::device_kokkos
{

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)

  /**
   * @brief Owning holder for one IPC mapping into another process's device allocation.
   *
   * Exists for two reasons that bare `void *` slots could not give:
   *  - every close runs through one place, so none of them can quietly discard its error code;
   *  - the compiler-generated move operations of the owner become correct. With bare pointers and
   *    defaulted moves, a moved-from ExchangeManager still held the pointers and its destructor
   *    closed mappings the moved-to object was using, while move-assignment leaked the target's.
   *
   * The slot is nulled *before* the close, so a slot never names an unmapped region -- not even
   * transiently, and not if the close fails.
   */
  class IpcMapping
  {
  public:
    IpcMapping() = default;
    IpcMapping(const IpcMapping &) = delete;
    IpcMapping &operator=(const IpcMapping &) = delete;
    IpcMapping(IpcMapping &&other) noexcept : mPtr(other.mPtr) { other.mPtr = nullptr; }
    IpcMapping &operator=(IpcMapping &&other) noexcept
    {
      if (this != &other) {
        reset();
        mPtr = other.mPtr;
        other.mPtr = nullptr;
      }
      return *this;
    }
    ~IpcMapping() { reset(); }

    void *get() const { return mPtr; }
    explicit operator bool() const { return mPtr != nullptr; }

    /** @brief Take ownership of a freshly opened mapping, closing whatever was held. */
    void adopt(void *ptr)
    {
      reset();
      mPtr = ptr;
    }

    /** @brief Close the mapping, if any, and null the slot. Never throws: this runs from
     *  destructors. A failing close is counted and reported instead of being discarded. */
    void reset() noexcept
    {
      if (mPtr == nullptr) return;
      void *ptr = mPtr;
      mPtr = nullptr;
      int err = p2p::ipcCloseHandle(ptr);
      if (err != 0) {
        ++p2p::ipcCloseErrorCount();
        try {
          sayMPI << "IPC close error: closing the mapping at " << ptr << " failed: " << p2p::ipcErrorString(err)
                 << ". The exporting rank most likely freed the allocation while it was still mapped here.\n";
        } catch (...) {
        }
      }
    }

  private:
    void *mPtr = nullptr;
  };

#endif

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
    ExchangeManager(MPICartesianExchange exchange, [[maybe_unused]] MPI_Comm shmComm, [[maybe_unused]] int myDeviceId)
        : mExchange(exchange)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      mMyDevice = myDeviceId;
      mCartComm = mExchange.getMPICartesianGroup().getComm();
      mShmComm = shmComm;
      MPI_Comm_rank(mCartComm, &mMyRank);

      mP2PAvailable.fill(false);
      mFullDuplex.fill(false);
      mRemoteHandleVersion.fill(0);

      // Neighbour ranks are recorded unconditionally, before any P2P probing: the handle
      // protocol below is collective over the whole cartesian communicator, so ranks with no
      // node-local peer of their own still take part and need their neighbours' ranks.
      auto &neighbours = mExchange.getNeighbours();
      for (size_t d = 0; d < NDim; ++d) {
        mNeighborRanks[d * 2 + 0] = neighbours.getUpperNeighbour(d);
        mNeighborRanks[d * 2 + 1] = neighbours.getLowerNeighbour(d);
      }

      probeP2P(shmComm);

      // Whether ANY rank in the job has a P2P link. The IPC handle protocol is entered on this
      // job-wide flag rather than on each rank's own links: its MPI_Sendrecv pairs are per-link
      // rendezvous, so a rank that opted out locally while a neighbour opted in would hang the
      // neighbour. This is the only predicate that is uniform by construction.
      int anyLocal = 0;
      for (size_t i = 0; i < 2 * NDim; ++i)
        anyLocal |= mP2PAvailable[i] ? 1 : 0;
      MPI_Allreduce(&anyLocal, &mAnyP2PInJob, 1, MPI_INT, MPI_MAX, mCartComm);
#endif
    }

    // Non-copyable (owns IPC mappings); movable -- IpcMapping nulls the source on move, so the
    // defaulted operations no longer double-close or leak.
    ExchangeManager(const ExchangeManager &) = delete;
    ExchangeManager &operator=(const ExchangeManager &) = delete;
    ExchangeManager(ExchangeManager &&) = default;
    ExchangeManager &operator=(ExchangeManager &&) = default;

    // ------------------------------------------------------------------
    // Buffer handle exchange — call after (re)allocating send/recv buffers
    // ------------------------------------------------------------------

    void updateBufferHandles([[maybe_unused]] char *sendUpPtr, [[maybe_unused]] char *sendDownPtr,
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
     * For P2P neighbors: barrier → P2P read from remote send buffer → streamSync → barrier.
     * For non-P2P neighbors: standard MPI Isend/Irecv (overlapped with P2P).
     */
    void exchange(size_t dimension, void *sendUpPtr, void *sendDownPtr, void *recvUpPtr, void *recvDownPtr,
                  size_t byteCount, int count, MPI_Datatype dataType)
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      // upP2P/downP2P describe the *neighbor*: upP2P = P2P-able with upper neighbor.
      // Sends and link-quality checks are gated by the neighbor flags. Receives use the
      // pull model — filling recvUp reads from the LOWER neighbor's sendUp, and filling
      // recvDown reads from the UPPER neighbor's sendDown — so they are gated by the
      // OPPOSITE neighbor's P2P-ability.
      bool upP2P = isP2PUp(dimension);
      bool downP2P = isP2PDown(dimension);
      bool canPullRecvUp = downP2P;
      bool canPullRecvDown = upP2P;

      if (upP2P || downP2P) {
        // Kokkos::fence was already called by the ghost updater (packing complete on GPU).
        bool allFullDuplex = (!upP2P || mFullDuplex[dimension * 2 + 0]) && (!downP2P || mFullDuplex[dimension * 2 + 1]);

        // Post non-P2P MPI receives before the handshake (overlaps with the token wait)
        if (!canPullRecvUp) mExchange.IrecvUp(dataType, dimension, recvUpPtr, count);
        if (!canPullRecvDown) mExchange.IrecvDown(dataType, dimension, recvDownPtr, count);

        // Pack-done handshake (replaces the pre-read shared-memory barrier): a pairwise 0-byte exchange
        // with each P2P neighbour of this dimension. It tells my readers my send buffers are packed, and
        // confirms the neighbours I pull from have packed theirs — the only ordering the old global
        // barrier actually provided for this rank. Device ordering still comes from the pack fence.
        p2pHandshake(dimension, upP2P, downP2P, MPITags::ghostP2PPackToken);

        // Post non-P2P MPI sends (can overlap with P2P reads)
        if (!upP2P) mExchange.IsendUp(dataType, dimension, sendUpPtr, count);
        if (!downP2P) mExchange.IsendDown(dataType, dimension, sendDownPtr, count);

        if (allFullDuplex) {
          // --- Single-phase: NVLink/xGMI is full-duplex, no bidirectional contention ---
          if (canPullRecvUp) p2p::memcpyAsync(recvUpPtr, mRemoteSendUpPtr[dimension * 2 + 1].get(), byteCount);
          if (canPullRecvDown) p2p::memcpyAsync(recvDownPtr, mRemoteSendDownPtr[dimension * 2 + 0].get(), byteCount);
          p2p::streamSynchronize();
        } else {
          // --- Two-phase: PCIe bidirectional contention avoidance ---
          // Simultaneous bidirectional P2P reads on a shared PCIe switch degrade throughput
          // by 10x+. Split reads by rank ordering: phase 0 if myRank < sourceRank,
          // phase 1 if myRank > sourceRank. No bidirectional pair in either phase.
          int upReadSource = mNeighborRanks[dimension * 2 + 1];   // lower neighbor
          int downReadSource = mNeighborRanks[dimension * 2 + 0]; // upper neighbor

          // Phase 0: reads where this rank has the lower rank number
          if (canPullRecvUp && mMyRank < upReadSource)
            p2p::memcpyAsync(recvUpPtr, mRemoteSendUpPtr[dimension * 2 + 1].get(), byteCount);
          if (canPullRecvDown && mMyRank < downReadSource)
            p2p::memcpyAsync(recvDownPtr, mRemoteSendDownPtr[dimension * 2 + 0].get(), byteCount);
          p2p::streamSynchronize();

          // Phase-ordering handshake (replaces the mid shared-memory barrier): on each P2P link the
          // lower-ranked rank reads in phase 0 then signals the higher-ranked rank, which waits for that
          // token before its phase-1 read — so the two reads on a link never overlap on the PCIe switch.
          p2pPhaseHandshake(dimension, canPullRecvUp, canPullRecvDown, upReadSource, downReadSource);

          // Phase 1: reads where this rank has the higher rank number
          if (canPullRecvUp && mMyRank > upReadSource)
            p2p::memcpyAsync(recvUpPtr, mRemoteSendUpPtr[dimension * 2 + 1].get(), byteCount);
          if (canPullRecvDown && mMyRank > downReadSource)
            p2p::memcpyAsync(recvDownPtr, mRemoteSendDownPtr[dimension * 2 + 0].get(), byteCount);
          p2p::streamSynchronize();
        }

        // Wait for non-P2P MPI to complete
        if (!upP2P || !downP2P) mExchange.waitall();

        // Read-done handshake (replaces the post-read shared-memory barrier): my send buffers are safe to
        // repack only once my P2P readers have finished pulling from them. Same pairwise neighbour set.
        p2pHandshake(dimension, upP2P, downP2P, MPITags::ghostP2PReadToken);

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

    /**
     * @brief Non-blocking in-place exchange of both up and down faces for one dimension (CPU path).
     *
     * Replaces the two sequential blocking MPI_Sendrecv (exchangeUp + exchangeDown) with four
     * concurrent Isend/Irecv and a single waitall, so the up and down transfers overlap.
     * Safe with the single dataShiftGhostCells tag: MPI non-overtaking matches (source,tag) in
     * posting order, and every rank posts recvs-before-sends and Up-before-Down consistently, so
     * even when a dimension has only two ranks (upper == lower neighbour) the matches are correct.
     * Send and receive regions never overlap (owned face vs ghost slice), so in-place is safe.
     * The dimension sweep stays sequential in GhostUpdater to preserve corner correctness.
     */
    void exchangeUpDownNonBlocking(MPI_Datatype dataType, ptrdiff_t dimension, void *sendUpPtr, void *recvUpPtr,
                                   void *sendDownPtr, void *recvDownPtr, int sendCount = 1)
    {
      mExchange.IrecvUp(dataType, dimension, recvUpPtr, sendCount);
      mExchange.IrecvDown(dataType, dimension, recvDownPtr, sendCount);
      mExchange.IsendUp(dataType, dimension, sendUpPtr, sendCount);
      mExchange.IsendDown(dataType, dimension, sendDownPtr, sendCount);
      mExchange.waitall();
    }

    /** @brief Coalesced non-blocking up/down exchange: each of the four messages carries its own
     *  absolute-address datatype and is addressed with MPI_BOTTOM. Used by the CPU component-coalesced
     *  ghost exchange to send all components of a multi-component field as one message per direction. */
    void exchangeUpDownBottom(ptrdiff_t dimension, MPI_Datatype sendUpType, MPI_Datatype recvUpType,
                              MPI_Datatype sendDownType, MPI_Datatype recvDownType)
    {
      mExchange.exchangeUpDownBottom(dimension, sendUpType, recvUpType, sendDownType, recvDownType);
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
    std::array<bool, 2 * NDim> mFullDuplex{}; // true if link is NVLink/xGMI (no bidirectional contention)
    // IPC-mapped pointers to each neighbor's send buffers (we READ from these)
    // mRemoteSendUpPtr[d*2+dir]: the neighbor in direction 'dir' of dimension d's sendUp buffer
    std::array<IpcMapping, 2 * NDim> mRemoteSendUpPtr{};
    std::array<IpcMapping, 2 * NDim> mRemoteSendDownPtr{};
    std::array<uint64_t, 2 * NDim> mRemoteHandleVersion{};
    std::array<int, 2 * NDim> mNeighborRanks{};
    std::array<int, 2 * NDim> mNeighborDevices{};
    int mAnyP2PInJob = 0; // job-wide: does any rank have a node-local P2P link?

    bool isP2PUp(size_t d) const { return mP2PAvailable[d * 2 + 0]; }
    bool isP2PDown(size_t d) const { return mP2PAvailable[d * 2 + 1]; }

    /**
     * @brief Pairwise 0-byte handshake with the (up to two) P2P neighbours of a dimension.
     *
     * Posts a non-blocking recv+send to each P2P neighbour, then a single Waitall -- deadlock-free for a
     * periodic dimension of ANY size. The earlier two ordered blocking MPI_Sendrecv (upper first, then
     * lower, on every rank) only avoided deadlock for size <= 2: on a ring of >= 3 same-node P2P ranks
     * every rank blocked in its upper Sendrecv waiting for the upper neighbour to reach its lower Sendrecv,
     * which never happened -> cyclic wait. Posting all recvs/sends before waiting removes that ordering
     * dependence. Token semantics are unchanged (mutual arrival on each link = my readers know my buffers
     * are packed / done being read). Relies on the P2P classification being symmetric (canAccessPeer is
     * symmetric in practice), which is also what the original barrier-based scheme required.
     */
    void p2pHandshake(size_t dimension, bool withUpper, bool withLower, int tag)
    {
      std::array<MPI_Request, 4> reqs;
      std::array<char, 2> sbuf{}, rbuf{}; // distinct buffers so the concurrent Irecvs never alias
      int n = 0;
      if (withUpper) {
        int up = mNeighborRanks[dimension * 2 + 0];
        MPI_Irecv(&rbuf[0], 1, MPI_BYTE, up, tag, mCartComm, &reqs[n++]);
        MPI_Isend(&sbuf[0], 1, MPI_BYTE, up, tag, mCartComm, &reqs[n++]);
      }
      if (withLower) {
        int lo = mNeighborRanks[dimension * 2 + 1];
        MPI_Irecv(&rbuf[1], 1, MPI_BYTE, lo, tag, mCartComm, &reqs[n++]);
        MPI_Isend(&sbuf[1], 1, MPI_BYTE, lo, tag, mCartComm, &reqs[n++]);
      }
      MPI_Waitall(n, reqs.data(), MPI_STATUSES_IGNORE);
    }

    /**
     * @brief Phase-ordering handshake for the two-phase (PCIe) P2P read, replacing the mid barrier.
     *
     * On each P2P link the lower-ranked rank reads in phase 0 and the higher-ranked rank in phase 1. After
     * its phase-0 read this rank sends a token to every source it out-ranks-below (mMyRank < source), and
     * waits for a token from every source it out-ranks-above (mMyRank > source) before starting phase 1.
     * So on any link the phase-1 read cannot begin until the phase-0 read on that same link has finished —
     * exactly the contention avoidance the barrier gave, but only between the two ranks that share the link.
     * Separate send/recv buffers keep concurrent Isend/Irecv (distinct sources) non-aliasing.
     */
    void p2pPhaseHandshake(size_t dimension, bool haveUp, bool haveDown, int upSource, int downSource)
    {
      std::array<MPI_Request, 4> reqs;
      std::array<char, 4> buf{}; // one distinct byte per in-flight message so concurrent recvs never alias
      int n = 0;
      const int tag = MPITags::ghostP2PPhaseToken;
      if (haveUp && mMyRank < upSource) {
        MPI_Isend(&buf[n], 1, MPI_BYTE, upSource, tag, mCartComm, &reqs[n]);
        ++n;
      }
      if (haveDown && mMyRank < downSource) {
        MPI_Isend(&buf[n], 1, MPI_BYTE, downSource, tag, mCartComm, &reqs[n]);
        ++n;
      }
      if (haveUp && mMyRank > upSource) {
        MPI_Irecv(&buf[n], 1, MPI_BYTE, upSource, tag, mCartComm, &reqs[n]);
        ++n;
      }
      if (haveDown && mMyRank > downSource) {
        MPI_Irecv(&buf[n], 1, MPI_BYTE, downSource, tag, mCartComm, &reqs[n]);
        ++n;
      }
      if (n > 0) MPI_Waitall(n, reqs.data(), MPI_STATUSES_IGNORE);
    }

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

      // mNeighborRanks was filled by the constructor, before this probe: it is needed whether or
      // not P2P turns out to be available.
      for (size_t d = 0; d < NDim; ++d) {
        checkAndEnableP2P(d, 0, mNeighborRanks[d * 2 + 0], rankDeviceMap);
        checkAndEnableP2P(d, 1, mNeighborRanks[d * 2 + 1], rankDeviceMap);
      }

      MPI_Group_free(&worldGroup);
      MPI_Group_free(&shmGroup);
    }

    void checkAndEnableP2P(size_t dim, int dir, int neighborRank, const std::vector<std::pair<int, int>> &rankDeviceMap)
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

      if (neighborDevice != mMyDevice) {
        if (!p2p::canAccessPeer(mMyDevice, neighborDevice)) return;
        p2p::enablePeerAccess(neighborDevice);
        mFullDuplex[idx] = p2p::isFullDuplexLink(mMyDevice, neighborDevice);
      } else {
        // Two ranks sharing the same GPU (e.g. GPU_NOCONSTRAIN oversubscription):
        // IPC between processes on the same device is valid and has no bus contention.
        // enablePeerAccess would error on self; skip it.
        mFullDuplex[idx] = true;
      }
      mP2PAvailable[idx] = true;

      sayMPI << "Ghost exchange: P2P enabled for dimension " << dim << (dir == 0 ? " (up)" : " (down)") << " to rank "
             << neighborRank << " (device " << neighborDevice << ", "
             << (neighborDevice == mMyDevice ? "same GPU" : (mFullDuplex[idx] ? "NVLink/xGMI" : "PCIe")) << ")\n";
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

      if (!mAnyP2PInJob) return;

      // The two handles do not depend on the dimension, so pack them once.
      p2p::IpcHandlePacket mySendUpPacket{};
      if (sendUpPtr != nullptr) p2p::ipcGetHandle(sendUpPtr, mySendUpPacket.handle);
      mySendUpPacket.deviceId = mMyDevice;
      mySendUpPacket.version = (sendUpPtr != nullptr) ? version : 0;

      p2p::IpcHandlePacket mySendDownPacket{};
      if (sendDownPtr != nullptr) p2p::ipcGetHandle(sendDownPtr, mySendDownPacket.handle);
      mySendDownPacket.deviceId = mMyDevice;
      mySendDownPacket.version = (sendDownPtr != nullptr) ? version : 0;

      // Discriminator, kept permanently because it costs one memcmp per growth and it names a
      // failure this design cannot survive: if the two send buffers are suballocated from a
      // single backing region the driver hands back the SAME handle for both, and an importer
      // that maps per (dimension, direction) opens it twice -> "resource already mapped". No
      // ordering of closes and opens fixes that; only deduplicating mappings per (peer,
      // allocation) does. Report it rather than abort, so a run that is otherwise healthy still
      // gets to say what went wrong.
      if (sendUpPtr != nullptr && sendDownPtr != nullptr &&
          std::memcmp(mySendUpPacket.handle, mySendDownPacket.handle, p2p::IpcHandleSize) == 0)
        sayMPI << "IPC handle alias: sendUp (" << static_cast<void *>(sendUpPtr) << ") and sendDown ("
               << static_cast<void *>(sendDownPtr) << ") produced identical IPC handles at version " << version
               << ". Importers keyed per (dimension, direction) will fail to open the second one.\n";

      for (size_t d = 0; d < NDim; ++d) {
        int upperRank = mNeighborRanks[d * 2 + 0];
        int lowerRank = mNeighborRanks[d * 2 + 1];

        // No gate on this rank's own P2P flags. The exchanges below are per-link rendezvous, so
        // a rank that skipped them while a neighbour with a node-local peer still posted its
        // half would hang that neighbour. The only skip every rank in the dimension agrees on
        // is an undivided dimension, where both neighbours are this rank itself.
        if (upperRank == mMyRank && lowerRank == mMyRank) continue;

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
        MPI_Sendrecv(&mySendDownPacket, sizeof(p2p::IpcHandlePacket), MPI_BYTE, lowerRank, tag2, &recvSendDownFromUpper,
                     sizeof(p2p::IpcHandlePacket), MPI_BYTE, upperRank, tag2, mCartComm, &stat);

        // Open lower neighbor's sendUp handle (for our recvUp)
        bool failLower = tryOpen(mRemoteSendUpPtr[d * 2 + 1], mP2PAvailable[d * 2 + 1], recvSendUpFromLower,
                                 mRemoteHandleVersion[d * 2 + 1], lowerRank, "sendUp");

        // Open upper neighbor's sendDown handle (for our recvDown)
        bool failUpper = tryOpen(mRemoteSendDownPtr[d * 2 + 0], mP2PAvailable[d * 2 + 0], recvSendDownFromUpper,
                                 mRemoteHandleVersion[d * 2 + 0], upperRank, "sendDown");

        // A failed open has to demote the whole LINK, on both ends. The same flag decides
        // whether the exporter posts an MPI send and whether the importer pulls, so a one-sided
        // fallback would leave one rank reading a mapping it does not have while the other
        // never sends. Each rank therefore tells a neighbour whether ITS open failed and learns
        // whether the neighbour's did -- same (dest, source) pairing as the handle exchange, so
        // what arrives from the lower neighbour concerns the link recorded in slot d*2+1.
        char myFailToUpper = failUpper ? 1 : 0;
        char myFailToLower = failLower ? 1 : 0;
        char peerFailFromLower = 0, peerFailFromUpper = 0;
        int tag3 = 700 + d * 4 + 2;
        MPI_Sendrecv(&myFailToUpper, 1, MPI_BYTE, upperRank, tag3, &peerFailFromLower, 1, MPI_BYTE, lowerRank, tag3,
                     mCartComm, &stat);
        int tag4 = 700 + d * 4 + 3;
        MPI_Sendrecv(&myFailToLower, 1, MPI_BYTE, lowerRank, tag4, &peerFailFromUpper, 1, MPI_BYTE, upperRank, tag4,
                     mCartComm, &stat);

        if (failLower || peerFailFromLower) demoteLink(d * 2 + 1, lowerRank);
        if (failUpper || peerFailFromUpper) demoteLink(d * 2 + 0, upperRank);

        mRemoteHandleVersion[d * 2 + 0] = std::max(mRemoteHandleVersion[d * 2 + 0], recvSendDownFromUpper.version);
        mRemoteHandleVersion[d * 2 + 1] = std::max(mRemoteHandleVersion[d * 2 + 1], recvSendUpFromLower.version);
      }
    }

    /**
     * @brief Open one published handle into `slot`. Returns true if the open FAILED.
     *
     * The slot is left empty on failure rather than left naming the region it named before --
     * the old code closed first and assigned second, so a throwing open left the slot pointing
     * at an unmapped region and the next exchange segfaulted inside cudaMemcpyAsync.
     */
    bool tryOpen(IpcMapping &slot, bool available, const p2p::IpcHandlePacket &packet, uint64_t knownVersion,
                 int peerRank, const char *what)
    {
      if (!available || packet.version == 0 || packet.version <= knownVersion) return false;
      slot.reset();
      try {
        slot.adopt(p2p::ipcOpenHandle(packet.handle));
      } catch (const std::exception &e) {
        sayMPI << "IPC open failed for rank " << peerRank << "'s " << what << " buffer (version " << packet.version
               << "): " << e.what() << ". Falling back to MPI on this link.\n";
        return true;
      }
      return false;
    }

    /** @brief Turn one P2P link off for the rest of the run and drop its mappings. Both ends do
     *  this in the same call, so the two ranks stay in agreement about who sends and who pulls. */
    void demoteLink(size_t idx, int peerRank)
    {
      mRemoteSendUpPtr[idx].reset();
      mRemoteSendDownPtr[idx].reset();
      if (!mP2PAvailable[idx]) return;
      mP2PAvailable[idx] = false;
      sayMPI << "Ghost exchange: P2P disabled for the link to rank " << peerRank
             << " after a failed IPC open; this link now uses MPI.\n";
    }
#endif
  };

} // namespace TempLat::device_kokkos

#endif // HAVE_MPI

#endif
