#ifndef TEMPLAT_LATTICE_MANIPULATION_GHOSTUPDATER_H
#define TEMPLAT_LATTICE_MANIPULATION_GHOSTUPDATER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler,  Year: 2025

#include "TempLat/parallel/mpi/mpitypeconstants.h"
#include "TempLat/parallel/mpi/mpitags.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesianexchange.h"
#include "TempLat/lattice/memory/memoryblock.h"
#include "TempLat/lattice/ghostcells/ghostsubarraymap.h"

#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"

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
        : mExchange(exchange), mLayout(layout), mGhostDepth(mLayout.getNGhosts()),
          mGhostSubarrayMap(mLayout, mGhostDepth)
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
          if (i != d)
            slabTotal *= full_sizes[i] + 2 * mGhostDepth;
        mMaxSlabSize = std::max(mMaxSlabSize, slabTotal);
      }
    }

    template <typename T> void update(MemoryBlock<T, NDim> &block)
    {
#ifdef HAVE_MPI
      // There is no MPI splitting in one dimension. Also, when we have only a single node, there is no need to do MPI
      // communication.
      if constexpr (NDim > 1) {
        if (mExchange.getMPICartesianGroup().size() > 1) {
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

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    MPICartesianExchange mExchange;
    LayoutStruct<NDim> mLayout;
    device::Idx mGhostDepth;
    GhostSubarrayMap<NDim> mGhostSubarrayMap;

    // Pre-computed max slab element count; buffers allocated lazily on first update<T>() call
    size_t mMaxSlabSize = 0;
    size_t mAllocatedBytes = 0;
    device::memory::NDView<char, 1> mSendUpBuffer;
    device::memory::NDView<char, 1> mRecvUpBuffer;
    device::memory::NDView<char, 1> mSendDownBuffer;
    device::memory::NDView<char, 1> mRecvDownBuffer;

    template <typename T> void pUpdate(MemoryBlock<T, NDim> &block)
    {
      /* iterate dimensions */
      for (size_t d = 0; d < NDim; ++d) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
        update_forDimension_device(block, d);
#else
        update_forDimension(block, d);
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
        mSendUpBuffer = device::memory::NDView<char, 1>("ghostSendUpBuf", neededBytes);
        mRecvUpBuffer = device::memory::NDView<char, 1>("ghostRecvUpBuf", neededBytes);
        mSendDownBuffer = device::memory::NDView<char, 1>("ghostSendDownBuf", neededBytes);
        mRecvDownBuffer = device::memory::NDView<char, 1>("ghostRecvDownBuf", neededBytes);
        mAllocatedBytes = neededBytes;
      }

      // Create unmanaged ND views over pre-allocated byte buffers (no cudaMalloc per call)
      auto sendUpSlab = device::apply(
          [&](const auto &...args) {
            return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(mSendUpBuffer.data()), args...);
          },
          slab_sizes);
      auto recvUpSlab = device::apply(
          [&](const auto &...args) {
            return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(mRecvUpBuffer.data()), args...);
          },
          slab_sizes);
      auto sendDownSlab = device::apply(
          [&](const auto &...args) {
            return device::memory::NDViewUnmanaged<T, NDim>(reinterpret_cast<T *>(mSendDownBuffer.data()), args...);
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
      auto sendUpSubView = device::apply(
          [&](const auto &...args) { return device::memory::subview(fullView, args...); }, sendUp_slices);
      auto sendDownSubView = device::apply(
          [&](const auto &...args) { return device::memory::subview(fullView, args...); }, sendDown_slices);
      device::memory::copyDeviceToDevice(sendUpSubView, sendUpSlab);
      device::memory::copyDeviceToDevice(sendDownSubView, sendDownSlab);
      device::iteration::fence(); // single fence ensures both packs complete

      // Launch all 4 non-blocking MPI operations
      MPI_Datatype dataType = MPITypeSelect<T>();
      mExchange.IrecvUp(dataType, dimension, recvUpSlab.data(), total_size);
      mExchange.IrecvDown(dataType, dimension, recvDownSlab.data(), total_size);
      mExchange.IsendUp(dataType, dimension, sendUpSlab.data(), total_size);
      mExchange.IsendDown(dataType, dimension, sendDownSlab.data(), total_size);
      mExchange.waitall();

      // Unpack both receive slabs (GPU kernels can run concurrently)
      auto recvUpSubView = device::apply(
          [&](const auto &...args) { return device::memory::subview(fullView, args...); }, recvUp_slices);
      auto recvDownSubView = device::apply(
          [&](const auto &...args) { return device::memory::subview(fullView, args...); }, recvDown_slices);
      device::memory::copyDeviceToDevice(recvUpSlab, recvUpSubView);
      device::memory::copyDeviceToDevice(recvDownSlab, recvDownSubView);
      device::iteration::fence(); // ensures both unpacks complete before next dimension
    }

  private:
    template <typename T> void update_forDimension(MemoryBlock<T, NDim> &block, device::Idx dimension)
    {
      auto *ptr = block.data();
#ifdef HAVE_MPI
      mExchange.exchangeUp(mGhostSubarrayMap.template getSubArray<T>(dimension), dimension,
                           /* base ptr is lower corner of all memory, including ghosts. */
                           /* send:
                            Don't jump to origin, but jump along the edge of dimension
                            to the point where we still have mGhostDepth until the end of
                            our *owned* memory (before the mGhostDepth hyper slices start) */
                           ptr + (mLayout.getSizesInMemory()[dimension]) * mLayout.stride(dimension),
                           /* receive: in origin, including ghosts. */
                           ptr);
      /* pointers: the same as above, but shifted by ghostDepth and ordering swapped. Yes. */
      mExchange.exchangeDown(mGhostSubarrayMap.template getSubArray<T>(dimension), dimension,
                             ptr + mGhostDepth * mLayout.stride(dimension),
                             ptr + (mGhostDepth + mLayout.getSizesInMemory()[dimension]) * mLayout.stride(dimension));
#endif
    }

  public:
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
