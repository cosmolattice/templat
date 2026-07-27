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
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
      device::p2p::rawDeviceFree(mSendUpRaw);
      device::p2p::rawDeviceFree(mSendDownRaw);
#else
      delete[] mSendUpRaw;
      delete[] mSendDownRaw;
#endif
    }

    template <typename T> void update(MemoryBlock<T, NDim> &block)
    {
      if (mGhostDepth == 0)
        throw GhostUpdaterException("Cannot update ghost cells with ghost depth 0. "
                                    "Use nGhost >= 1 when creating MemoryToolBox.");
#ifdef HAVE_MPI
      // There is no MPI splitting in one dimension. Also, when we have only a single node, there is no need to do MPI
      // communication.
      if constexpr (NDim > 1) {
        if (mExchangeManager.getMPICartesianGroup().size() > 1) {
          pUpdate(block);
        } else {
          pUpdate_NOMPI(block);
        }
      } else
#endif
      {
        pUpdate_NOMPI(block);
      }
    }

    /** @brief Coalesced ghost update of several equally-shaped blocks (the components of a
     *  multi-component field). All faces of a given dimension/direction are packed into one buffer and
     *  exchanged in a single message, cutting the per-component message and synchronization count.
     *  The dimension sweep stays sequential so corner/edge ghosts remain correct. A batch of one is
     *  byte-identical to update(). */
    template <typename T> void updateBatch(std::span<MemoryBlock<T, NDim> *const> blocks)
    {
      if (blocks.empty()) return;
      if (blocks.size() == 1) {
        update(*blocks[0]);
        return;
      }
      if (mGhostDepth == 0)
        throw GhostUpdaterException("Cannot update ghost cells with ghost depth 0. "
                                    "Use nGhost >= 1 when creating MemoryToolBox.");
#ifdef HAVE_MPI
      if constexpr (NDim > 1) {
        if (mExchangeManager.getMPICartesianGroup().size() > 1) {
          pUpdateBatch(blocks);
        } else {
          for (auto *b : blocks)
            pUpdate_NOMPI(*b);
        }
      } else
#endif
      {
        for (auto *b : blocks)
          pUpdate_NOMPI(*b);
      }
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
#ifdef HAVE_MPI
    device::memory::ExchangeManager<NDim> mExchangeManager;
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

    template <typename T> void pUpdate(MemoryBlock<T, NDim> &block)
    {
#ifdef HAVE_MPI
      auto &decomp = mExchangeManager.getMPICartesianGroup().getDecomposition();
#endif
      /* iterate dimensions */
      for (size_t d = 0; d < NDim; ++d) {
#ifdef HAVE_MPI
        // Non-split dimensions: local periodic copy (no MPI overhead)
        if (decomp[d] <= 1) {
          pUpdate_NOMPI_singleDim(block, d);
          continue;
        }
#endif
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
        update_forDimension_device(block, d);
#else
        update_forDimension(block, d);
#endif
      }
    }

    /** @brief Coalesced multi-block variant of pUpdate. Keeps the dimension sweep outermost and
     *  sequential (corner correctness); within a split dimension all blocks' faces are exchanged
     *  together (one MPI message / one P2P pull per direction) on both the CPU and GPU paths. */
    template <typename T> void pUpdateBatch(std::span<MemoryBlock<T, NDim> *const> blocks)
    {
#ifdef HAVE_MPI
      auto &decomp = mExchangeManager.getMPICartesianGroup().getDecomposition();
#endif
      for (size_t d = 0; d < NDim; ++d) {
#ifdef HAVE_MPI
        // Non-split dimensions: local periodic copy per block (no MPI overhead).
        if (decomp[d] <= 1) {
          for (auto *b : blocks)
            pUpdate_NOMPI_singleDim(*b, d);
          continue;
        }
#endif
#if defined(DEVICE_CUDA) || defined(DEVICE_HIP)
        update_forDimension_device_batch(blocks, d);
#else
        update_forDimension_batch(blocks, d);
#endif
      }
    }

  public:
    template <typename T> void update_forDimension_device(MemoryBlock<T, NDim> &block, size_t dimension)
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
    }

    /** @brief Coalesced GPU exchange of one split dimension for several blocks. Each component's face is
     *  packed into the shared slab buffers at its own contiguous offset (c * total_size), so the whole
     *  batch travels as one exchange (one P2P pull / one MPI message of C*total_size elements) with a
     *  single pack fence and a single unpack fence, instead of C separate exchanges + 2C fences. The
     *  send buffers grow to C * maxSlab and re-publish their IPC handles once per new maximum C. */
    template <typename T>
    void update_forDimension_device_batch(std::span<MemoryBlock<T, NDim> *const> blocks, size_t dimension)
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
    }

  private:
    template <typename T> void update_forDimension(MemoryBlock<T, NDim> &block, device::Idx dimension)
    {
      auto *ptr = block.data();
#ifdef HAVE_MPI
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
    void update_forDimension_batch(std::span<MemoryBlock<T, NDim> *const> blocks, device::Idx dimension)
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
#endif
    }

  public:
    /** @brief Local periodic ghost copy for a single dimension (no MPI). */
    template <typename T> void pUpdate_NOMPI_singleDim(MemoryBlock<T, NDim> &block, size_t dim)
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
        if constexpr (NDim == 1) {
          device::iteration::foreach (
              "GhostUpdater", device::IdxArray<1>{0}, device::IdxArray<1>{1},
              DEVICE_LAMBDA(const device::IdxArray<1> &i) {
                View(ghostDepth - depth) = View(ghostDepth + sizes[0] - depth);
                View(ghostDepth + sizes[0] + (depth - 1)) = View(ghostDepth + (depth - 1));
              });
        } else {
          device::array<std::pair<device::Idx, device::Idx>, NDim> btf_slicesFrom{};
          device::array<std::pair<device::Idx, device::Idx>, NDim> btf_slicesTo{};
          device::array<std::pair<device::Idx, device::Idx>, NDim> ftb_slicesFrom{};
          device::array<std::pair<device::Idx, device::Idx>, NDim> ftb_slicesTo{};

          for (size_t i = 0; i < NDim; ++i) {
            btf_slicesFrom[i] = (i == dim)
                                    ? std::make_pair<device::Idx, device::Idx>(ghostDepth + sizes[i] - depth,
                                                                               ghostDepth + sizes[i] - depth + 1)
                                    : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
            btf_slicesTo[i] = (i == dim)
                                  ? std::make_pair<device::Idx, device::Idx>(ghostDepth - depth, ghostDepth - depth + 1)
                                  : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
            ftb_slicesFrom[i] =
                (i == dim)
                    ? std::make_pair<device::Idx, device::Idx>(ghostDepth + (depth - 1), ghostDepth + (depth - 1) + 1)
                    : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
            ftb_slicesTo[i] = (i == dim)
                                  ? std::make_pair<device::Idx, device::Idx>(ghostDepth + sizes[i] + (depth - 1),
                                                                             ghostDepth + sizes[i] + (depth - 1) + 1)
                                  : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
          }
          auto btf_fromSubView = device::apply(
              [&](const auto &...args) { return device::memory::subview(View, args...); }, btf_slicesFrom);
          auto btf_toSubView =
              device::apply([&](const auto &...args) { return device::memory::subview(View, args...); }, btf_slicesTo);
          auto ftb_fromSubView = device::apply(
              [&](const auto &...args) { return device::memory::subview(View, args...); }, ftb_slicesFrom);
          auto ftb_toSubView =
              device::apply([&](const auto &...args) { return device::memory::subview(View, args...); }, ftb_slicesTo);

          device::memory::copyDeviceToDevice(btf_fromSubView, btf_toSubView);
          device::memory::copyDeviceToDevice(ftb_fromSubView, ftb_toSubView);
        }
      }
    }

    template <typename T> void pUpdate_NOMPI(MemoryBlock<T, NDim> &block, device::Idx dimension = 0)
    {
      // Get View to the full data
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

      // Create subviews for the from and to views
      // We need to create slices for each dimension, taking into account the padding
      // and the layout of the views
      device::array<std::pair<device::Idx, device::Idx>, NDim> btf_slicesFrom{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> btf_slicesTo{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> ftb_slicesFrom{};
      device::array<std::pair<device::Idx, device::Idx>, NDim> ftb_slicesTo{};

      for (size_t dim = 0; dim < NDim; ++dim) {
        for (size_t depth = 1; depth <= (size_t)mGhostDepth; ++depth) {

          if constexpr (NDim == 1) {
            // For NDim == 1, we just need to copy the corners.
            device::iteration::foreach (
                "GhostUpdater", device::IdxArray<1>{0}, device::IdxArray<1>{1},
                DEVICE_LAMBDA(const device::IdxArray<1> &i) {
                  View(ghostDepth - depth) = View(ghostDepth + sizes[0] - depth);
                  View(ghostDepth + sizes[0] + (depth - 1)) = View(ghostDepth + (depth - 1));
                });
          } else {
            // so we copy a (NDim- 1)-dimensional slice. Include the padding, which leads to a copy of all corners, too!
            for (size_t i = 0; i < NDim; ++i) {
              btf_slicesFrom[i] = (i == dim)
                                      ? std::make_pair<device::Idx, device::Idx>(ghostDepth + sizes[i] - depth,
                                                                                 ghostDepth + sizes[i] - depth + 1)
                                      : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
              btf_slicesTo[i] =
                  (i == dim) ? std::make_pair<device::Idx, device::Idx>(ghostDepth - depth, ghostDepth - depth + 1)
                             : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
              ftb_slicesFrom[i] =
                  (i == dim)
                      ? std::make_pair<device::Idx, device::Idx>(ghostDepth + (depth - 1), ghostDepth + (depth - 1) + 1)
                      : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
              ftb_slicesTo[i] = (i == dim)
                                    ? std::make_pair<device::Idx, device::Idx>(ghostDepth + sizes[i] + (depth - 1),
                                                                               ghostDepth + sizes[i] + (depth - 1) + 1)
                                    : std::make_pair<device::Idx, device::Idx>(0, ghostDepth + sizes[i] + ghostDepth);
            }
            auto btf_fromSubView = device::apply(
                [&](const auto &...args) { return device::memory::subview(View, args...); }, btf_slicesFrom);
            auto btf_toSubView = device::apply(
                [&](const auto &...args) { return device::memory::subview(View, args...); }, btf_slicesTo);
            auto ftb_fromSubView = device::apply(
                [&](const auto &...args) { return device::memory::subview(View, args...); }, ftb_slicesFrom);
            auto ftb_toSubView = device::apply(
                [&](const auto &...args) { return device::memory::subview(View, args...); }, ftb_slicesTo);

            // Copy the data in both directions
            device::memory::copyDeviceToDevice(btf_fromSubView, btf_toSubView);
            device::memory::copyDeviceToDevice(ftb_fromSubView, ftb_toSubView);
          }
        }
      }
    }
  };
} // namespace TempLat

#endif
