#ifndef TEMPLAT_LATTICE_MANIPULATION_GHOSTUPDATER_H
#define TEMPLAT_LATTICE_MANIPULATION_GHOSTUPDATER_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2025

#include "TempLat/parallel/mpi/mpitypeconstants.h"
#include "TempLat/parallel/mpi/mpitags.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesianexchange.h"
#include "TempLat/lattice/memory/memoryblock.h"
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/ghostcells/ghostsubarraymap.h"

#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/parallel/device_guard.h"

#include <span>
#include <vector>

namespace TempLat
{
  MakeException(GhostUpdaterException);

  /** @brief A class which updates the ghost cells in our total memory block.
   * By having the LayoutStruct, this class knows what is the ghostDepth.
   *
   * Has one public method, update<T>(T* ptr), which, based on LayoutStruct,
   * uses the associated subarrays and performs the exchange up and down
   * in all dimensions, through calls to MPICartesianExchange with the
   * appriate datatypes for the subarrays.
   *
   *
   * Unit test: ctest -R test-ghostupdater
   **/
  template <size_t NDim> class GhostUpdater
  {
  public:
    // Put public methods here. These should change very little over time.
    GhostUpdater(MPICartesianExchange exchange, LayoutStruct<NDim> layout)
        :
#ifdef HAVE_MPI
          mExchangeManager(exchange, DeviceGuard::getShmComm(), DeviceGuard::getDeviceId()),
#endif
          mLayout(layout), mGhostDepth(mLayout.getNGhosts()), mGhostSubarrayMap(mLayout, mGhostDepth)
    {
      auto full_sizes = mLayout.getSizesInMemory();
      for (size_t i = 0; i < NDim; ++i) {
        if (mGhostDepth > full_sizes[i]) {
          throw GhostUpdaterException("Ghost depth is larger than local size in dimension " + std::to_string(i) + ":",
                                      mGhostDepth, " > ", full_sizes[i]);
        }
      }
      /* verify that */
      bool allSame = true;
      for (auto &&it : mLayout.getPadding()) {
        allSame = allSame && mGhostDepth == it[0] && mGhostDepth == it[1];
      }
      if (!allSame)
        throw GhostUpdaterException(
            "Can only work with identical padding at start and end of each dimension, not this.", allSame);

      // Pre-allocate GPU slab buffers for ghost exchange to avoid per-call cudaMalloc/cudaFree
      mMaxSlabSize = 0;
      for (size_t d = 0; d < NDim; ++d) {
        size_t slabTotal = mGhostDepth;
        for (size_t i = 0; i < NDim; ++i)
          if (i != d) slabTotal *= full_sizes[i] + 2 * mGhostDepth;
        mMaxSlabSize = std::max(mMaxSlabSize, slabTotal);
      }
    }

    ~GhostUpdater()
    {
#ifdef HAVE_MPI
      // Teardown is the other place these send buffers are freed, so it needs the same
      // close-before-free ordering as growth does; without it every rank frees buffers its
      // peers may still have mapped and the closes that follow report on freed memory. Safe as
      // a collective for the same reason destruction of the shared MemoryToolBox is: all ranks
      // reach it. A rank unwinding alone would already be heading for MPI_Abort.
      mExchangeManager.retireBufferHandles();
#endif
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
      device::p2p::rawDeviceFree(mSendUpRaw);
      device::p2p::rawDeviceFree(mSendDownRaw);
#else
      delete[] mSendUpRaw;
      delete[] mSendDownRaw;
#endif
    }

    template <typename T>
    void update(MemoryBlock<T, NDim> &block, BCSpec<NDim> bcSpec = allPeriodic<NDim>())
    {
      if (mGhostDepth == 0)
        throw GhostUpdaterException("Cannot update ghost cells with ghost depth 0. "
                                    "Use nGhost >= 1 when creating MemoryToolBox.");
#ifdef HAVE_MPI
      // There is no MPI splitting in one dimension. Also, when we have only a single node, there is no need to do MPI
      // communication.
      if constexpr (NDim > 1) {
        if (mExchangeManager.getMPICartesianGroup().size() > 1) {
          pUpdate(block, bcSpec);
        } else {
          pUpdate_NOMPI(block, bcSpec);
        }
      } else
#endif
      {
        pUpdate_NOMPI(block, bcSpec);
      }
    }

    /** @brief Coalesced ghost update of several equally-shaped blocks (the components of a
     *  multi-component field). All faces of a given dimension/direction are packed into one buffer and
     *  exchanged in a single message, cutting the per-component message and synchronization count.
     *  The dimension sweep stays sequential so corner/edge ghosts remain correct. A batch of one is
     *  byte-identical to update().
     *
     *  This overload applies one `bcSpec` to every block in the batch; see the per-component
     *  overload below for fields whose components differ. */
    template <typename T>
    void updateBatch(std::span<MemoryBlock<T, NDim> *const> blocks,
                     BCSpec<NDim> bcSpec = allPeriodic<NDim>())
    {
      if (blocks.empty()) return;
      std::vector<BCSpec<NDim>> bcSpecs(blocks.size(), bcSpec);
      updateBatch(blocks, std::span<const BCSpec<NDim>>(bcSpecs.data(), bcSpecs.size()));
    }

    /** @brief updateBatch with per-component boundary conditions: `bcSpecs[c]` belongs to
     *  `blocks[c]`.
     *
     *  Components of one field disagreeing on their BC is not pathological — it is what C-star
     *  boundary conditions are. The identification Phi(x + L) = Phi*(x) is, in the real-component
     *  basis these fields are stored in, a sign flip on some components and not others:
     *  (c0, c1, c2, c3) -> (c0, -c1, c2, -c3) for both an SU(2) group element and an SU(2)
     *  doublet, (Re, Im) -> (Re, -Im) for a complex scalar. So a batch legitimately carries a
     *  mixture of Periodic and Antiperiodic.
     *
     *  Nothing about the coalescing depends on uniformity: the exchange is BC-agnostic (it just
     *  lands the global-wrap value in every ghost slab) and the BC is a per-block post-step that
     *  was already written as a loop over the batch. Only the fused local kernel needs care —
     *  see pUpdate_NOMPI_singleDim_batch. */
    template <typename T>
    void updateBatch(std::span<MemoryBlock<T, NDim> *const> blocks,
                     std::span<const BCSpec<NDim>> bcSpecs)
    {
      if (blocks.empty()) return;
      if (bcSpecs.size() != blocks.size())
        throw GhostUpdaterException("updateBatch: expected exactly one BCSpec per block.");
      if (blocks.size() == 1) {
        update(*blocks[0], bcSpecs[0]);
        return;
      }
      if (mGhostDepth == 0)
        throw GhostUpdaterException("Cannot update ghost cells with ghost depth 0. "
                                    "Use nGhost >= 1 when creating MemoryToolBox.");
#ifdef HAVE_MPI
      if constexpr (NDim > 1) {
        if (mExchangeManager.getMPICartesianGroup().size() > 1) {
          pUpdateBatch(blocks, bcSpecs);
        } else {
          pUpdate_NOMPI_batchLocal(blocks, bcSpecs);
        }
      } else
#endif
      {
        pUpdate_NOMPI_batchLocal(blocks, bcSpecs);
      }
    }

  private:
    /** @brief The single-rank leg of updateBatch. Coalescing across components used to stop at
     *  the MPI boundary: with one rank (or with every dimension undecomposed) updateBatch fell
     *  through to a per-component pUpdate_NOMPI loop, so a 4-component field paid 4x the ghost
     *  launches even though the batch machinery was right there. On one GPU that was ~1174 copy
     *  launches per MC sweep against 70 sweep kernels. */
    template <typename T>
    void pUpdate_NOMPI_batchLocal(std::span<MemoryBlock<T, NDim> *const> blocks,
                                  std::span<const BCSpec<NDim>> bcSpecs)
    {
      pUpdate_NOMPI_batch<T>(blocks, bcSpecs);
    }

    /** @brief Does any component of the batch ask for a non-Periodic BC in `dim`? Lets the MPI
     *  boundary post-step keep its single all-or-nothing early-out for the (overwhelmingly
     *  common) all-periodic batch, now that the BC is read per component inside the loop. */
    static bool anyNonPeriodic(std::span<const BCSpec<NDim>> bcSpecs, size_t dim)
    {
      for (const auto &s : bcSpecs)
        if (s[dim] != BCType::Periodic) return true;
      return false;
    }

  public:

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
#ifdef HAVE_MPI
    device::memory::ExchangeManager<NDim> mExchangeManager;

    /** @brief Returns {isLowBoundary, isHighBoundary} for this rank along `dim`. */
    std::pair<bool, bool> isBoundaryRank(size_t dim) const
    {
      const auto &group = mExchangeManager.getMPICartesianGroup();
      const int coord = group.getPosition()[dim];
      const int decomp = group.getDecomposition()[dim];
      return {coord == 0, coord == decomp - 1};
    }
#endif
    LayoutStruct<NDim> mLayout;
    device::Idx mGhostDepth;
    GhostSubarrayMap<NDim> mGhostSubarrayMap;

    // Pre-computed max slab element count; buffers allocated lazily on first update<T>() call
    size_t mMaxSlabSize = 0;
    size_t mAllocatedBytes = 0;
    uint64_t mHandleVersion = 0;

    // Send buffers: raw GPU allocations so IPC handles point to exact data start.
    // Recv buffers: views (local-only, no IPC needed).
    char *mSendUpRaw = nullptr;
    char *mSendDownRaw = nullptr;
    device::memory::NDView<char, 1> mRecvUpBuffer;
    device::memory::NDView<char, 1> mRecvDownBuffer;

    template <typename T>
    void pUpdate(MemoryBlock<T, NDim> &block, BCSpec<NDim> bcSpec = allPeriodic<NDim>())
    {
#ifdef HAVE_MPI
      auto &decomp = mExchangeManager.getMPICartesianGroup().getDecomposition();
#endif
      /* iterate dimensions */
      for (size_t d = 0; d < NDim; ++d) {
#ifdef HAVE_MPI
        // Non-split dimensions: local periodic copy (no MPI overhead)
        if (decomp[d] <= 1) {
          pUpdate_NOMPI_singleDim(block, d, bcSpec);
          continue;
        }
#endif
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
        update_forDimension_device(block, d, bcSpec);
#else
        update_forDimension(block, d, bcSpec);
#endif
      }
    }

    /** @brief Coalesced multi-block variant of pUpdate. Keeps the dimension sweep outermost and
     *  sequential (corner correctness); within a split dimension all blocks' faces are exchanged
     *  together (one MPI message / one P2P pull per direction) on both the CPU and GPU paths. */
    template <typename T>
    void pUpdateBatch(std::span<MemoryBlock<T, NDim> *const> blocks,
                      std::span<const BCSpec<NDim>> bcSpecs)
    {
#ifdef HAVE_MPI
      auto &decomp = mExchangeManager.getMPICartesianGroup().getDecomposition();
#endif
      for (size_t d = 0; d < NDim; ++d) {
#ifdef HAVE_MPI
        // Non-split dimensions: local BC-aware copy, no MPI. Batched across components for the
        // same reason the split dimensions are -- with a 1-D decomposition two of three
        // dimensions land here, so this is most of the ghost work even under MPI.
        if (decomp[d] <= 1) {
          pUpdate_NOMPI_singleDim_batch<T>(blocks, d, bcSpecs);
          continue;
        }
#endif
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
        update_forDimension_device_batch(blocks, d, bcSpecs);
#else
        update_forDimension_batch(blocks, d, bcSpecs);
#endif
      }
    }

  public:
    template <typename T>
    void update_forDimension_device(MemoryBlock<T, NDim> &block, size_t dimension,
                                    BCSpec<NDim> bcSpec = allPeriodic<NDim>())
    {
      // We will copy slabs of thickness ghostDepth in the dimension 'dimension'.
      device::IdxArray<NDim> full_sizes = mLayout.getSizesInMemory();
      for (size_t i = 0; i < NDim; ++i)
        full_sizes[i] += 2 * mGhostDepth;
      device::IdxArray<NDim> slab_sizes = mLayout.getSizesInMemory();
      for (size_t i = 0; i < NDim; ++i)
        slab_sizes[i] += 2 * mGhostDepth;
      slab_sizes[dimension] = mGhostDepth;
      size_t total_size = 1;
      for (size_t i = 0; i < NDim; ++i)
        total_size *= slab_sizes[i];

      // Ensure byte buffers are large enough for this T (lazy alloc on first call or type change)
      size_t neededBytes = mMaxSlabSize * sizeof(T);
      if (neededBytes > mAllocatedBytes) {
#ifdef HAVE_MPI
        // Retire, barrier, reallocate, republish. Peers hold IPC mappings into these send
        // buffers, so every rank must close its mappings BEFORE any rank frees -- see
        // ExchangeManager::retireBufferHandles.
        //
        // This is a collective call, which is only safe because growth happens in lockstep:
        // mAllocatedBytes is always mMaxSlabSize times the largest sizeof(T)*C seen so far, and
        // mMaxSlabSize is a positive per-rank constant, so the predicate above reduces to
        // sizeof(T)*C > that largest product -- identical on every rank even when the
        // decomposition is uneven. All ranks therefore take this branch in the same call.
        mExchangeManager.retireBufferHandles();
#endif
        // Send buffers: raw GPU alloc for clean IPC base pointers
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
        device::p2p::rawDeviceFree(mSendUpRaw);
        device::p2p::rawDeviceFree(mSendDownRaw);
        mSendUpRaw = static_cast<char *>(device::p2p::rawDeviceMalloc(neededBytes));
        mSendDownRaw = static_cast<char *>(device::p2p::rawDeviceMalloc(neededBytes));
#else
        // CPU fallback: use operator new
        delete[] mSendUpRaw;
        delete[] mSendDownRaw;
        mSendUpRaw = new char[neededBytes];
        mSendDownRaw = new char[neededBytes];
#endif
        // Recv buffers: views (local-only)
        mRecvUpBuffer = device::memory::NDView<char, 1>("ghostRecvUpBuf", neededBytes);
        mRecvDownBuffer = device::memory::NDView<char, 1>("ghostRecvDownBuf", neededBytes);
        mAllocatedBytes = neededBytes;
#ifdef HAVE_MPI
        mExchangeManager.updateBufferHandles(mSendUpRaw, mSendDownRaw, ++mHandleVersion);
#endif
      }

      // Create unmanaged ND views over pre-allocated byte buffers (no cudaMalloc per call)
      auto sendUpSlab = device::apply(
          [&](const auto &...args) {
            return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(mSendUpRaw), args...);
          },
          slab_sizes);
      auto recvUpSlab = device::apply(
          [&](const auto &...args) {
            return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(mRecvUpBuffer.data()), args...);
          },
          slab_sizes);
      auto sendDownSlab = device::apply(
          [&](const auto &...args) {
            return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(mSendDownRaw), args...);
          },
          slab_sizes);
      auto recvDownSlab = device::apply(
          [&](const auto &...args) {
            return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(mRecvDownBuffer.data()), args...);
          },
          slab_sizes);

      // Compute slices for UP and DOWN directions
      device::array<std::pair<device::Idx, device::Idx>, NDim> sendUp_slices{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> recvUp_slices{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> sendDown_slices{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> recvDown_slices{};

      for (size_t i = 0; i < NDim; ++i) {
        // UP: send end of dimension, receive at origin
        sendUp_slices[i] = (i == dimension) ? std::pair<device::Idx, device::Idx>(full_sizes[i] - 2 * mGhostDepth,
                                                                                  full_sizes[i] - mGhostDepth)
                                            : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
        recvUp_slices[i] = (i == dimension) ? std::pair<device::Idx, device::Idx>(0, mGhostDepth)
                                            : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
        // DOWN: send origin of dimension, receive at end
        sendDown_slices[i] = (i == dimension) ? std::pair<device::Idx, device::Idx>(mGhostDepth, 2 * mGhostDepth)
                                              : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
        recvDown_slices[i] = (i == dimension)
                                 ? std::pair<device::Idx, device::Idx>(full_sizes[i] - mGhostDepth, full_sizes[i])
                                 : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
      }

      auto fullView = block.getNDView(full_sizes);

      // Pack both UP and DOWN send slabs (GPU kernels can run concurrently)
      auto sendUpSubView =
          device::apply([&](const auto &...args) { return device::memory::subview(fullView, args...); }, sendUp_slices);
      auto sendDownSubView = device::apply(
          [&](const auto &...args) { return device::memory::subview(fullView, args...); }, sendDown_slices);
      device::memory::copyDeviceToDevice(sendUpSubView, sendUpSlab);
      device::memory::copyDeviceToDevice(sendDownSubView, sendDownSlab);
      device::iteration::fence(); // single fence ensures both packs complete

      // Exchange ghost slabs — ExchangeManager routes to P2P or MPI per direction
#ifdef HAVE_MPI
      MPI_Datatype dataType = MPITypeSelect<T>();
      mExchangeManager.exchange(dimension, sendUpSlab.data(), sendDownSlab.data(), recvUpSlab.data(),
                                recvDownSlab.data(), total_size * sizeof(T), total_size, dataType);
#endif

      // Unpack both receive slabs (GPU kernels can run concurrently)
      auto recvUpSubView =
          device::apply([&](const auto &...args) { return device::memory::subview(fullView, args...); }, recvUp_slices);
      auto recvDownSubView = device::apply(
          [&](const auto &...args) { return device::memory::subview(fullView, args...); }, recvDown_slices);
      device::memory::copyDeviceToDevice(recvUpSlab, recvUpSubView);
      device::memory::copyDeviceToDevice(recvDownSlab, recvDownSubView);
      device::iteration::fence(); // ensures both unpacks complete before next dimension

#ifdef HAVE_MPI
      // BC-aware post-step: on boundary ranks along `dimension`, overwrite the wrap-around data
      // just unpacked into the low-/high-ghost slab with the BC transform. Non-boundary ranks and
      // Periodic BC are untouched (the exchange already produced correct values).
      // recvUpSubView is the LOW-ghost slab (slice [0, mGhostDepth) in dim) — filled from lower
      //                neighbor; at coord==0 this is wrap-around from the far end of the Cart.
      // recvDownSubView is the HIGH-ghost slab (slice [full-mGhostDepth, full) in dim) — filled
      //                from upper neighbor; at coord==decomp-1 this is wrap-around.
      if (bcSpec[dimension] != BCType::Periodic) {
        const auto boundary = isBoundaryRank(dimension);
        if (boundary.first || boundary.second) {
          device::IdxArray<NDim> ownedSizes = mLayout.getSizesInMemory();
          for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth) {
            applyLocalBCAtDimDepth<T>(fullView, dimension, depth, ownedSizes, mGhostDepth,
                                      bcSpec[dimension], boundary.first, boundary.second,
                                      /*mpiPostStep=*/true);
          }
        }
      }
#else
      (void)bcSpec;
#endif
    }

    /** @brief Coalesced GPU exchange of one split dimension for several blocks. Each component's face is
     *  packed into the shared slab buffers at its own contiguous offset (c * total_size), so the whole
     *  batch travels as one exchange (one P2P pull / one MPI message of C*total_size elements) with a
     *  single pack fence and a single unpack fence, instead of C separate exchanges + 2C fences. The
     *  send buffers grow to C * maxSlab and re-publish their IPC handles once per new maximum C. */
    template <typename T>
    void update_forDimension_device_batch(std::span<MemoryBlock<T, NDim> *const> blocks, size_t dimension,
                                          std::span<const BCSpec<NDim>> bcSpecs)
    {
      const size_t C = blocks.size();

      device::IdxArray<NDim> full_sizes = mLayout.getSizesInMemory();
      for (size_t i = 0; i < NDim; ++i)
        full_sizes[i] += 2 * mGhostDepth;
      device::IdxArray<NDim> slab_sizes = mLayout.getSizesInMemory();
      for (size_t i = 0; i < NDim; ++i)
        slab_sizes[i] += 2 * mGhostDepth;
      slab_sizes[dimension] = mGhostDepth;
      size_t total_size = 1;
      for (size_t i = 0; i < NDim; ++i)
        total_size *= slab_sizes[i];

      // Grow the shared send/recv buffers to hold all C component slabs contiguously.
      size_t neededBytes = mMaxSlabSize * sizeof(T) * C;
      if (neededBytes > mAllocatedBytes) {
#ifdef HAVE_MPI
        // Collective close-before-free; see the identical comment in update_forDimension_device
        // for why this branch is guaranteed to be taken by every rank in the same call.
        mExchangeManager.retireBufferHandles();
#endif
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
        device::p2p::rawDeviceFree(mSendUpRaw);
        device::p2p::rawDeviceFree(mSendDownRaw);
        mSendUpRaw = static_cast<char *>(device::p2p::rawDeviceMalloc(neededBytes));
        mSendDownRaw = static_cast<char *>(device::p2p::rawDeviceMalloc(neededBytes));
#else
        delete[] mSendUpRaw;
        delete[] mSendDownRaw;
        mSendUpRaw = new char[neededBytes];
        mSendDownRaw = new char[neededBytes];
#endif
        mRecvUpBuffer = device::memory::NDView<char, 1>("ghostRecvUpBuf", neededBytes);
        mRecvDownBuffer = device::memory::NDView<char, 1>("ghostRecvDownBuf", neededBytes);
        mAllocatedBytes = neededBytes;
#ifdef HAVE_MPI
        mExchangeManager.updateBufferHandles(mSendUpRaw, mSendDownRaw, ++mHandleVersion);
#endif
      }

      // Slices are identical for every component (same layout).
      device::array<std::pair<device::Idx, device::Idx>, NDim> sendUp_slices{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> recvUp_slices{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> sendDown_slices{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> recvDown_slices{};
      for (size_t i = 0; i < NDim; ++i) {
        sendUp_slices[i] = (i == dimension) ? std::pair<device::Idx, device::Idx>(full_sizes[i] - 2 * mGhostDepth,
                                                                                  full_sizes[i] - mGhostDepth)
                                            : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
        recvUp_slices[i] = (i == dimension) ? std::pair<device::Idx, device::Idx>(0, mGhostDepth)
                                            : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
        sendDown_slices[i] = (i == dimension) ? std::pair<device::Idx, device::Idx>(mGhostDepth, 2 * mGhostDepth)
                                              : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
        recvDown_slices[i] = (i == dimension)
                                 ? std::pair<device::Idx, device::Idx>(full_sizes[i] - mGhostDepth, full_sizes[i])
                                 : std::pair<device::Idx, device::Idx>(0, slab_sizes[i]);
      }

      // Unmanaged slab view over one component's slot inside a byte buffer.
      auto slabView = [&](char *base, size_t c) {
        return device::apply(
            [&](const auto &...args) {
              return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(base) + c * total_size, args...);
            },
            slab_sizes);
      };

      // Pack every component's up/down faces into its slot. One fence for the whole batch.
      for (size_t c = 0; c < C; ++c) {
        auto fullView = blocks[c]->getNDView(full_sizes);
        auto sendUpSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(fullView, args...); }, sendUp_slices);
        auto sendDownSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(fullView, args...); }, sendDown_slices);
        auto sendUpSlab = slabView(mSendUpRaw, c);
        auto sendDownSlab = slabView(mSendDownRaw, c);
        device::memory::copyDeviceToDevice(sendUpSubView, sendUpSlab);
        device::memory::copyDeviceToDevice(sendDownSubView, sendDownSlab);
      }
      device::iteration::fence();

      // One coalesced exchange of all C slabs (C*total_size elements).
#ifdef HAVE_MPI
      MPI_Datatype dataType = MPITypeSelect<T>();
      mExchangeManager.exchange(dimension, mSendUpRaw, mSendDownRaw, mRecvUpBuffer.data(), mRecvDownBuffer.data(),
                                C * total_size * sizeof(T), static_cast<int>(C * total_size), dataType);
#endif

      // Unpack every component from its slot. One fence for the whole batch.
      for (size_t c = 0; c < C; ++c) {
        auto fullView = blocks[c]->getNDView(full_sizes);
        auto recvUpSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(fullView, args...); }, recvUp_slices);
        auto recvDownSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(fullView, args...); }, recvDown_slices);
        device::memory::copyDeviceToDevice(slabView(mRecvUpBuffer.data(), c), recvUpSubView);
        device::memory::copyDeviceToDevice(slabView(mRecvDownBuffer.data(), c), recvDownSubView);
      }
      device::iteration::fence();

#ifdef HAVE_MPI
      // BC-aware post-step, per component — mirrors update_forDimension_device()'s fixup, looped over
      // the batch. Placed after the unpack fence so every component's ghost slab holds the exchanged
      // wrap-around value before the BC transform rewrites it in place. The BC is read per component:
      // the exchange above was BC-agnostic, so components carrying different BCs cost nothing extra
      // here beyond the ones that are Periodic being skipped.
      if (anyNonPeriodic(bcSpecs, dimension)) {
        const auto boundary = isBoundaryRank(dimension);
        if (boundary.first || boundary.second) {
          device::IdxArray<NDim> ownedSizes = mLayout.getSizesInMemory();
          for (size_t c = 0; c < C; ++c) {
            if (bcSpecs[c][dimension] == BCType::Periodic) continue;
            auto fullView = blocks[c]->getNDView(full_sizes);
            for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth) {
              applyLocalBCAtDimDepth<T>(fullView, dimension, depth, ownedSizes, mGhostDepth,
                                        bcSpecs[c][dimension], boundary.first, boundary.second,
                                        /*mpiPostStep=*/true);
            }
          }
        }
      }
#else
      (void)bcSpecs;
#endif
    }

  private:
    template <typename T>
    void update_forDimension(MemoryBlock<T, NDim> &block, device::Idx dimension,
                             BCSpec<NDim> bcSpec = allPeriodic<NDim>())
    {
#ifdef HAVE_MPI
      auto *ptr = block.data();
      // Non-blocking exchange of both faces of this dimension: the up and down transfers are posted
      // concurrently (Isend/Irecv) and awaited once, rather than as two sequential blocking MPI_Sendrecv.
      // The dimension sweep in pUpdate() stays sequential, so corner cells (which rely on previously
      // filled dimensions' ghosts being included) remain correct. Pointers match the previous
      // exchangeUp/exchangeDown layout:
      //   up   : send the top owned slice   -> recv into the lower ghost (origin)
      //   down : send the first owned slice -> recv into the upper ghost
      auto subArray = mGhostSubarrayMap.template getSubArray<T>(dimension);
      auto *sendUpPtr = ptr + (mLayout.getSizesInMemory()[dimension]) * mLayout.stride(dimension);
      auto *recvUpPtr = ptr;
      auto *sendDownPtr = ptr + mGhostDepth * mLayout.stride(dimension);
      auto *recvDownPtr = ptr + (mGhostDepth + mLayout.getSizesInMemory()[dimension]) * mLayout.stride(dimension);
      mExchangeManager.exchangeUpDownNonBlocking(subArray, dimension, sendUpPtr, recvUpPtr, sendDownPtr, recvDownPtr);

      // BC-aware post-step: on boundary ranks, overwrite the wrap-around ghost slab with the BC
      // transform. Uses the same applyLocalBCAtDimDepth helper as Phase 2 / device path. Periodic
      // BC and non-boundary ranks are untouched. exchangeUpDownNonBlocking waitall()s internally,
      // so the exchanged data is already in place here.
      if (bcSpec[dimension] != BCType::Periodic) {
        const auto boundary = isBoundaryRank(dimension);
        if (boundary.first || boundary.second) {
          device::IdxArray<NDim> ownedSizes = mLayout.getSizesInMemory();
          device::IdxArray<NDim> full_sizes{};
          for (size_t i = 0; i < NDim; ++i) full_sizes[i] = ownedSizes[i] + 2 * mGhostDepth;
          auto fullView = block.getNDView(full_sizes);
          for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth) {
            applyLocalBCAtDimDepth<T>(fullView, dimension, depth, ownedSizes, mGhostDepth,
                                      bcSpec[dimension], boundary.first, boundary.second,
                                      /*mpiPostStep=*/true);
          }
        }
      }
#else
      (void)block;
      (void)dimension;
      (void)bcSpec;
#endif
    }

    /** @brief Coalesced CPU exchange of one split dimension for several blocks. Instead of C separate
     *  strided-subarray messages per direction, we build one MPI datatype per direction that gathers all
     *  C component faces at their absolute addresses (MPI_Type_create_hindexed_block over the cached
     *  per-dimension subarray type), and send it with MPI_BOTTOM. MPI performs exactly the same strided
     *  gather/scatter it did per component, but as a single up-message and single down-message — so the
     *  transport does no extra copy and the per-component message/latency count collapses to one. The
     *  dimension sweep stays sequential (corner correctness). */
    template <typename T>
    void update_forDimension_batch(std::span<MemoryBlock<T, NDim> *const> blocks, device::Idx dimension,
                                   std::span<const BCSpec<NDim>> bcSpecs)
    {
#ifdef HAVE_MPI
      auto subArrayHolder = mGhostSubarrayMap.template getSubArray<T>(dimension);
      MPI_Datatype subArray = subArrayHolder;

      const size_t C = blocks.size();
      const auto stride = mLayout.stride(dimension);
      const auto owned = mLayout.getSizesInMemory()[dimension];

      // Absolute byte addresses of each component's face slab, per direction. Pointers match the single
      // block layout: up sends the top owned slice and recvs into the lower ghost (origin); down sends
      // the first owned slice and recvs into the upper ghost.
      std::vector<MPI_Aint> dispSendUp(C), dispRecvUp(C), dispSendDown(C), dispRecvDown(C);
      for (size_t c = 0; c < C; ++c) {
        T *ptr = blocks[c]->data();
        MPI_Get_address(ptr + owned * stride, &dispSendUp[c]);
        MPI_Get_address(ptr, &dispRecvUp[c]);
        MPI_Get_address(ptr + mGhostDepth * stride, &dispSendDown[c]);
        MPI_Get_address(ptr + (mGhostDepth + owned) * stride, &dispRecvDown[c]);
      }

      const int n = static_cast<int>(C);
      MPI_Datatype tSendUp, tRecvUp, tSendDown, tRecvDown;
      MPI_Type_create_hindexed_block(n, 1, dispSendUp.data(), subArray, &tSendUp);
      MPI_Type_create_hindexed_block(n, 1, dispRecvUp.data(), subArray, &tRecvUp);
      MPI_Type_create_hindexed_block(n, 1, dispSendDown.data(), subArray, &tSendDown);
      MPI_Type_create_hindexed_block(n, 1, dispRecvDown.data(), subArray, &tRecvDown);
      MPI_Type_commit(&tSendUp);
      MPI_Type_commit(&tRecvUp);
      MPI_Type_commit(&tSendDown);
      MPI_Type_commit(&tRecvDown);

      // One up-message and one down-message, both faces concurrent, buffers addressed via MPI_BOTTOM.
      mExchangeManager.exchangeUpDownBottom(dimension, tSendUp, tRecvUp, tSendDown, tRecvDown);

      MPI_Type_free(&tSendUp);
      MPI_Type_free(&tRecvUp);
      MPI_Type_free(&tSendDown);
      MPI_Type_free(&tRecvDown);

      // BC-aware post-step, per block. The coalesced exchange above is BC-agnostic: it lands the
      // global-wrap value in every component's ghost slab exactly as the single-block path does, so
      // the fixup is the same one update_forDimension() applies, just looped over the batch. Runs
      // inside the per-dimension call so the sequential dimension sweep (corner correctness) holds.
      // Read per component so a C-star batch (a mix of Periodic and Antiperiodic components) is
      // handled by the same loop.
      if (anyNonPeriodic(bcSpecs, dimension)) {
        const auto boundary = isBoundaryRank(dimension);
        if (boundary.first || boundary.second) {
          device::IdxArray<NDim> ownedSizes = mLayout.getSizesInMemory();
          device::IdxArray<NDim> full_sizes{};
          for (size_t i = 0; i < NDim; ++i) full_sizes[i] = ownedSizes[i] + 2 * mGhostDepth;
          for (size_t c = 0; c < C; ++c) {
            if (bcSpecs[c][dimension] == BCType::Periodic) continue;
            auto fullView = blocks[c]->getNDView(full_sizes);
            for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth) {
              applyLocalBCAtDimDepth<T>(fullView, dimension, depth, ownedSizes, mGhostDepth,
                                        bcSpecs[c][dimension], boundary.first, boundary.second,
                                        /*mpiPostStep=*/true);
            }
          }
        }
      }
#else
      (void)blocks;
      (void)dimension;
      (void)bcSpecs;
#endif
    }

    // Public, not private, only because of an nvcc restriction: the function enclosing an
    // extended __host__ __device__ lambda (which is what DEVICE_LAMBDA expands to under CUDA)
    // must not have private or protected access within its class. Both helpers below are
    // implementation details of the BC path and are not meant to be called from outside.
  public:
    // Elementwise dst = -src over all sites. Same shape required. GPU-safe via foreach.
    // Guarded by if-constexpr on unary minus so the Antiperiodic branch of applyLocalBCAtDimDepth
    // stays instantiable for element types that do not support negation (e.g. test-only struct types
    // that only exercise the default Periodic path).
    template <typename SrcView, typename DstView>
    void negatingCopySubview(const SrcView &src, DstView &dst)
    {
      static_assert(DstView::rank == SrcView::rank, "negatingCopySubview: rank mismatch");
      using ValT = typename DstView::value_type;
      if constexpr (requires(ValT v) { -v; }) {
        constexpr size_t R = DstView::rank;
        device::IdxArray<R> ends;
        for (size_t i = 0; i < R; ++i) ends[i] = dst.extent(i);
        device::iteration::foreach (
            "GhostUpdaterAntiperiodic", device::IdxArray<R>{}, ends,
            DEVICE_LAMBDA(const device::IdxArray<R> &idx) {
              device::apply([&](const auto... i) { dst(i...) = -src(i...); }, idx);
            });
      } else {
        throw GhostUpdaterException(
            "Antiperiodic BC requires element type to support unary operator-.");
      }
    }

    // BC-aware local ghost fill at a single (dim, depth). Writes the low- and/or high-ghost slabs
    // of `view` along `dim` at the requested `depth` level according to `bc`:
    //   Periodic:    wrap from opposite face.
    //   Antiperiodic: wrap from opposite face with sign flip.
    //   Dirichlet:   zero-fill (no source read).
    //   Neumann:     mirror the innermost interior cells outward (even extension).
    // The `doLow` / `doHigh` flags select which side(s) to write. The local (single-rank) path
    // calls with both true; the MPI boundary-rank post-step calls with exactly one true.
    //
    // `mpiPostStep` distinguishes the local single-rank path (false) from the MPI boundary-rank
    // post-exchange step (true). Under MPI, the exchange has already populated the destination
    // ghost slab with the global-wrap value; for Antiperiodic the BC is then an in-place sign
    // flip of that slab, NOT a copy from the local last/first owned cell (which under MPI lives
    // on the wrong end of the global domain). Dirichlet/Neumann/Periodic produce identical
    // results either way on a boundary rank — Dirichlet doesn't read source, Neumann's local
    // first/last IS the global first/last on a boundary rank, and Periodic skips the post-step.
    template <typename T, typename View>
    void applyLocalBCAtDimDepth(View &view, size_t dim, size_t depth,
                                const device::IdxArray<NDim> &sizes, device::Idx ghostDepth, BCType bc,
                                bool doLow = true, bool doHigh = true, bool mpiPostStep = false)
    {
      if (!doLow && !doHigh) return;
      if constexpr (NDim == 1) {
        switch (bc) {
        case BCType::Periodic:
          device::iteration::foreach (
              "GhostUpdater", device::IdxArray<1>{0}, device::IdxArray<1>{1},
              DEVICE_LAMBDA(const device::IdxArray<1> &) {
                if (doLow)  view(ghostDepth - depth) = view(ghostDepth + sizes[0] - depth);
                if (doHigh) view(ghostDepth + sizes[0] + (depth - 1)) = view(ghostDepth + (depth - 1));
              });
          break;
        case BCType::Antiperiodic:
          if (mpiPostStep) {
            device::iteration::foreach (
                "GhostUpdater", device::IdxArray<1>{0}, device::IdxArray<1>{1},
                DEVICE_LAMBDA(const device::IdxArray<1> &) {
                  if (doLow)  view(ghostDepth - depth) = -view(ghostDepth - depth);
                  if (doHigh) view(ghostDepth + sizes[0] + (depth - 1)) = -view(ghostDepth + sizes[0] + (depth - 1));
                });
          } else {
            device::iteration::foreach (
                "GhostUpdater", device::IdxArray<1>{0}, device::IdxArray<1>{1},
                DEVICE_LAMBDA(const device::IdxArray<1> &) {
                  if (doLow)  view(ghostDepth - depth) = -view(ghostDepth + sizes[0] - depth);
                  if (doHigh) view(ghostDepth + sizes[0] + (depth - 1)) = -view(ghostDepth + (depth - 1));
                });
          }
          break;
        case BCType::Dirichlet:
          device::iteration::foreach (
              "GhostUpdater", device::IdxArray<1>{0}, device::IdxArray<1>{1},
              DEVICE_LAMBDA(const device::IdxArray<1> &) {
                if (doLow)  view(ghostDepth - depth) = T{0};
                if (doHigh) view(ghostDepth + sizes[0] + (depth - 1)) = T{0};
              });
          break;
        case BCType::Neumann:
          device::iteration::foreach (
              "GhostUpdater", device::IdxArray<1>{0}, device::IdxArray<1>{1},
              DEVICE_LAMBDA(const device::IdxArray<1> &) {
                if (doLow)  view(ghostDepth - depth) = view(ghostDepth + (depth - 1));
                if (doHigh) view(ghostDepth + sizes[0] + (depth - 1)) = view(ghostDepth + sizes[0] - depth);
              });
          break;
        }
        return;
      } else {
        // Destination slices are BC-independent: low ghost at (ghostDepth - depth), high ghost
        // at (ghostDepth + sizes[dim] + (depth - 1)). Source slices differ:
        //   Periodic/Antiperiodic: opposite face (wrap).
        //   Neumann:              inward mirror (swap roles vs. wrap).
        //   Dirichlet:            no source needed — fill destination with T{0}.
        const bool mirror = (bc == BCType::Neumann);

        device::array<std::pair<device::Idx, device::Idx>, NDim> btf_slicesFrom{};
        device::array<std::pair<device::Idx, device::Idx>, NDim> btf_slicesTo{};
        device::array<std::pair<device::Idx, device::Idx>, NDim> ftb_slicesFrom{};
        device::array<std::pair<device::Idx, device::Idx>, NDim> ftb_slicesTo{};

        for (size_t i = 0; i < NDim; ++i) {
          const auto fullOther =
              std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
          btf_slicesTo[i] =
              (i == dim)
                  ? std::make_pair<device::Idx, device::Idx>(ghostDepth - depth, ghostDepth - depth + 1)
                  : fullOther;
          ftb_slicesTo[i] =
              (i == dim)
                  ? std::make_pair<device::Idx, device::Idx>(ghostDepth + sizes[i] + (depth - 1),
                                                             ghostDepth + sizes[i] + (depth - 1) + 1)
                  : fullOther;
          btf_slicesFrom[i] =
              (i == dim)
                  ? (mirror ? std::make_pair<device::Idx, device::Idx>(ghostDepth + (depth - 1),
                                                                       ghostDepth + (depth - 1) + 1)
                            : std::make_pair<device::Idx, device::Idx>(ghostDepth + sizes[i] - depth,
                                                                       ghostDepth + sizes[i] - depth + 1))
                  : fullOther;
          ftb_slicesFrom[i] =
              (i == dim)
                  ? (mirror ? std::make_pair<device::Idx, device::Idx>(ghostDepth + sizes[i] - depth,
                                                                       ghostDepth + sizes[i] - depth + 1)
                            : std::make_pair<device::Idx, device::Idx>(ghostDepth + (depth - 1),
                                                                       ghostDepth + (depth - 1) + 1))
                  : fullOther;
        }

        auto btf_fromSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(view, args...); }, btf_slicesFrom);
        auto btf_toSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(view, args...); }, btf_slicesTo);
        auto ftb_fromSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(view, args...); }, ftb_slicesFrom);
        auto ftb_toSubView = device::apply(
            [&](const auto &...args) { return device::memory::subview(view, args...); }, ftb_slicesTo);

        switch (bc) {
        case BCType::Periodic:
        case BCType::Neumann:
          if (doLow)  device::memory::copyDeviceToDevice(btf_fromSubView, btf_toSubView);
          if (doHigh) device::memory::copyDeviceToDevice(ftb_fromSubView, ftb_toSubView);
          break;
        case BCType::Antiperiodic:
          if (mpiPostStep) {
            // Sign-flip the exchanged ghost slab in place: ghost = -ghost. The exchange has
            // already brought the global-wrap value into the destination; we just negate it.
            if (doLow)  negatingCopySubview(btf_toSubView, btf_toSubView);
            if (doHigh) negatingCopySubview(ftb_toSubView, ftb_toSubView);
          } else {
            if (doLow)  negatingCopySubview(btf_fromSubView, btf_toSubView);
            if (doHigh) negatingCopySubview(ftb_fromSubView, ftb_toSubView);
          }
          break;
        case BCType::Dirichlet:
          if (doLow)  device::memory::fill(btf_toSubView, T{0});
          if (doHigh) device::memory::fill(ftb_toSubView, T{0});
          break;
        }
      }
    }

    /** @brief Whether `T` supports unary minus, i.e. whether antiperiodic BCs are expressible
     *  for it at all. Only the fused kernel needs to ask: the per-component path reaches its
     *  negation through a branch that a non-negatable type never instantiates. */
    template <typename U>
    static constexpr bool cCanNegate = requires(const U &x) { -x; };

    /** @brief `neg ? -v : v`, with the negation instantiated only when `U` can be negated.
     *
     *  A free function rather than a `if constexpr` inside the kernel body: nvcc rejects an
     *  extended __host__ __device__ lambda that *first-captures* a variable inside a
     *  constexpr-if context, which is exactly what capturing `negate` there amounts to
     *  ("cannot first-capture variable in constexpr-if context"). Confining the branch to a
     *  named DEVICE_INLINE_FUNCTION keeps the lambda's captures unconditional.
     *
     *  The `else` branch is unreachable in practice -- applyLocalBCAtDimDepthBatch throws
     *  before launch if antiperiodic BCs are asked of a non-negatable type -- so it exists to
     *  make the kernel compile for such types, not to define behaviour for them. */
    template <typename U>
    static DEVICE_INLINE_FUNCTION U negateIf(const U &v, bool neg)
    {
      if constexpr (cCanNegate<U>) return neg ? -v : v;
      else return v;
    }

    /** @brief Largest component count the fused local ghost kernel handles in one launch.
     *  The widest batching caller today is SymTracelessField (5 components); SU2Field and
     *  SU2Doublet use 4. Wider batches fall back to the per-component path, so this is a
     *  capture-size budget rather than a correctness limit. */
    static constexpr size_t cMaxLocalGhostBatch = 8;

    /** @brief Fused local BC ghost fill at one (dim, depth) across a batch of components.
     *
     *  applyLocalBCAtDimDepth issues two strided copyDeviceToDevice per component per
     *  (dim, depth) -- a low-face and a high-face copy. Over a full MC sweep that is ~1174
     *  copies, against 70 actual sweep kernels, and at small volumes the launch count rather
     *  than the byte count is what costs. This does the same writes for every component in
     *  ONE kernel: 6C copies per pass collapse to 3 launches (one per dimension).
     *
     *  Unlike the MPI batch path there is no slab to pack -- a local wrap is a direct
     *  device-to-device move -- so the saving has to come from kernel fusion, not from
     *  coalescing messages. Indices are therefore computed in-kernel instead of via subviews,
     *  which also keeps the captured argument small (one view array, not four subview arrays).
     *
     *  Only the non-mpiPostStep semantics are implemented, which is all the local path needs;
     *  the MPI boundary post-step keeps using the per-component routine.
     */
    template <typename T>
    void applyLocalBCAtDimDepthBatch(std::span<MemoryBlock<T, NDim> *const> blocks, size_t dim,
                                     size_t depth, const device::IdxArray<NDim> &sizes,
                                     device::Idx ghostDepth, BCType bc)
    {
      const size_t nb = blocks.size();

      device::IdxArray<NDim> full_sizes{};
      for (size_t i = 0; i < NDim; ++i)
        full_sizes[i] = ghostDepth + sizes[i] + ghostDepth;

      using ViewT = decltype(blocks[0]->template getNDView<T>(full_sizes));
      device::array<ViewT, cMaxLocalGhostBatch> views{};
      for (size_t c = 0; c < nb; ++c)
        views[c] = blocks[c]->template getNDView<T>(full_sizes);

      // Destination planes are BC-independent; only the source differs. Mirrors exactly the
      // btf_/ftb_ slice arithmetic in applyLocalBCAtDimDepth.
      const bool mirror = (bc == BCType::Neumann);
      const device::Idx dstLow = ghostDepth - depth;
      const device::Idx dstHigh = ghostDepth + sizes[dim] + (depth - 1);
      const device::Idx srcLow = mirror ? ghostDepth + (depth - 1) : ghostDepth + sizes[dim] - depth;
      const device::Idx srcHigh = mirror ? ghostDepth + sizes[dim] - depth : ghostDepth + (depth - 1);
      const bool negate = (bc == BCType::Antiperiodic);
      const bool zero = (bc == BCType::Dirichlet);

      // Antiperiodic is the only BC that needs unary minus, and the per-component path only
      // ever instantiates it inside its Antiperiodic branch (negatingCopySubview). The fused
      // kernel below applies it through negateIf(), which is where the requirement is confined;
      // written inline as `negate ? -v : v` it would instead demand operator- from EVERY
      // component type under EVERY boundary condition -- not a requirement the path it replaces
      // imposes, and one that does not hold for types carrying no arithmetic at all
      // (test-ghostupdaterbatch's labelled<N> markers being the case that caught it).
      // Refuse here the one combination negateIf() cannot honour, rather than silently
      // dropping a sign.
      if constexpr (!cCanNegate<T>)
        if (negate)
          throw GhostUpdaterException("Antiperiodic boundary conditions require a component type "
                                      "supporting unary minus; this one does not.");

      // Iterate the face: full extent (ghosts included, which is what fills corners) in every
      // direction but `dim`, collapsed to a single plane along `dim`. Sources and destinations
      // are disjoint -- interior vs ghost -- so writing both faces in one kernel is safe.
      device::array<device::Idx, NDim> starts{};
      device::array<device::Idx, NDim> stops{};
      for (size_t i = 0; i < NDim; ++i) {
        starts[i] = 0;
        stops[i] = (i == dim) ? 1 : full_sizes[i];
      }

      device::iteration::foreach (
          "GhostUpdaterBatch", starts, stops,
          DEVICE_LAMBDA(const device::IdxArray<NDim> &idx) {
            device::IdxArray<NDim> lo = idx, hi = idx, sl = idx, sh = idx;
            lo[dim] = dstLow;
            hi[dim] = dstHigh;
            sl[dim] = srcLow;
            sh[dim] = srcHigh;
            for (size_t c = 0; c < nb; ++c) {
              if (zero) {
                device::apply([&](const auto &...a) { views[c](a...) = T{0}; }, lo);
                device::apply([&](const auto &...a) { views[c](a...) = T{0}; }, hi);
              } else {
                const auto vlo = device::apply([&](const auto &...a) { return views[c](a...); }, sl);
                const auto vhi = device::apply([&](const auto &...a) { return views[c](a...); }, sh);
                device::apply([&](const auto &...a) { views[c](a...) = negateIf(vlo, negate); }, lo);
                device::apply([&](const auto &...a) { views[c](a...) = negateIf(vhi, negate); }, hi);
              }
            }
          });
    }

    /** @brief Fused local ghost update of one dimension across a batch of components. */
    template <typename T>
    void pUpdate_NOMPI_singleDim_batch(std::span<MemoryBlock<T, NDim> *const> blocks, size_t dim,
                                       std::span<const BCSpec<NDim>> bcSpecs)
    {
      if (blocks.empty()) return;
      // NDim == 1 has its own hand-written branch in applyLocalBCAtDimDepth, and an
      // over-wide batch would blow the capture budget: defer both to the per-component path.
      if constexpr (NDim == 1) {
        for (size_t c = 0; c < blocks.size(); ++c)
          pUpdate_NOMPI_singleDim(*blocks[c], dim, bcSpecs[c]);
        return;
      } else {
        if (blocks.size() > cMaxLocalGhostBatch) {
          for (size_t c = 0; c < blocks.size(); ++c)
            pUpdate_NOMPI_singleDim(*blocks[c], dim, bcSpecs[c]);
          return;
        }
        const auto ghostDepth = mLayout.getPadding()[0][0];
        device::IdxArray<NDim> sizes;
        for (size_t i = 0; i < NDim; ++i)
          sizes[i] = mLayout.getSizesInMemory()[i];

        // The fused kernel below writes one BCType for the whole batch: the source plane it reads
        // and the transform it applies are both BC-dependent, so they cannot be per-component
        // without pushing a per-component branch into the inner loop. Instead, partition the batch
        // by the BC it carries *in this dimension* and launch once per distinct BCType present.
        // A uniform batch (every field except a C-star one) takes the fast path below unchanged
        // and is byte-identical to before; a C-star batch has two groups, so 2 launches instead of
        // 1 -- still far below the C launches the per-component path would cost.
        bool uniform = true;
        for (size_t c = 1; c < bcSpecs.size(); ++c)
          if (bcSpecs[c][dim] != bcSpecs[0][dim]) { uniform = false; break; }

        if (uniform) {
          for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth)
            applyLocalBCAtDimDepthBatch<T>(blocks, dim, depth, sizes, ghostDepth, bcSpecs[0][dim]);
          return;
        }

        for (BCType bc : {BCType::Periodic, BCType::Antiperiodic, BCType::Dirichlet, BCType::Neumann}) {
          std::vector<MemoryBlock<T, NDim> *> group;
          for (size_t c = 0; c < blocks.size(); ++c)
            if (bcSpecs[c][dim] == bc) group.push_back(blocks[c]);
          if (group.empty()) continue;
          auto groupSpan = std::span<MemoryBlock<T, NDim> *const>(group.data(), group.size());
          for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth)
            applyLocalBCAtDimDepthBatch<T>(groupSpan, dim, depth, sizes, ghostDepth, bc);
        }
      }
    }

    /** @brief Fused local ghost update over all dimensions for a batch of components.
     *
     *  The dimension sweep stays outermost and sequential, which is what makes corner and edge
     *  ghosts correct: dimension d reads the ghosts that dimension d-1 just wrote. Note this
     *  inverts the loop nesting relative to `for (b : blocks) pUpdate_NOMPI(*b)`, which is
     *  component-outer / dimension-inner. That is safe because components are independent
     *  allocations -- what must be ordered is the dimensions *within* a component, and every
     *  component still sees them in the same order. */
    template <typename T>
    void pUpdate_NOMPI_batch(std::span<MemoryBlock<T, NDim> *const> blocks,
                             std::span<const BCSpec<NDim>> bcSpecs)
    {
      for (size_t d = 0; d < NDim; ++d)
        pUpdate_NOMPI_singleDim_batch<T>(blocks, d, bcSpecs);
    }

  public:
    /** @brief Local BC-aware ghost copy for a single dimension (no MPI). */
    template <typename T>
    void pUpdate_NOMPI_singleDim(MemoryBlock<T, NDim> &block, size_t dim,
                                 BCSpec<NDim> bcSpec = allPeriodic<NDim>())
    {
      const auto ghostDepth = mLayout.getPadding()[0][0];
      device::IdxArray<NDim> sizes;
      for (size_t i = 0; i < NDim; ++i)
        sizes[i] = mLayout.getSizesInMemory()[i];
      device::IdxArray<NDim> full_sizes{};
      for (size_t i = 0; i < NDim; ++i)
        full_sizes[i] = ghostDepth + sizes[i] + ghostDepth;
      auto View = block.getNDView(full_sizes);

      for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth) {
        applyLocalBCAtDimDepth<T>(View, dim, depth, sizes, ghostDepth, bcSpec[dim]);
      }
    }

    template <typename T>
    void pUpdate_NOMPI(MemoryBlock<T, NDim> &block, BCSpec<NDim> bcSpec = allPeriodic<NDim>(),
                      device::Idx dimension = 0)
    {
      const auto ghostDepth = mLayout.getPadding()[0][0];
      for (size_t i = 0; i < NDim; ++i)
        if (ghostDepth != mLayout.getPadding()[i][0] || ghostDepth != mLayout.getPadding()[i][1]) {
          throw GhostUpdaterException(
              "Can only work with identical padding at start and end of each dimension, not this.");
        }
      device::IdxArray<NDim> sizes;
      for (size_t i = 0; i < NDim; ++i)
        sizes[i] = mLayout.getSizesInMemory()[i];
      device::IdxArray<NDim> full_sizes{};
      for (size_t i = 0; i < NDim; ++i)
        full_sizes[i] = ghostDepth + sizes[i] + ghostDepth;
      auto View = block.getNDView(full_sizes);

      for (size_t dim = 0; dim < NDim; ++dim) {
        for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth) {
          applyLocalBCAtDimDepth<T>(View, dim, depth, sizes, ghostDepth, bcSpec[dim]);
        }
      }
    }
  };
} // namespace TempLat

#endif
