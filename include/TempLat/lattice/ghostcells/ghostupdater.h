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
#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/lattice/ghostcells/ghostsubarraymap.h"

#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/parallel/device_guard.h"

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

  private:
    template <typename T>
    void update_forDimension(MemoryBlock<T, NDim> &block, device::Idx dimension,
                             BCSpec<NDim> bcSpec = allPeriodic<NDim>())
    {
#ifdef HAVE_MPI
      auto *ptr = block.data();
      mExchangeManager.exchangeUp(mGhostSubarrayMap.template getSubArray<T>(dimension), dimension,
                                  /* base ptr is lower corner of all memory, including ghosts. */
                                  /* send:
                                   Don't jump to origin, but jump along the edge of dimension
                                   to the point where we still have mGhostDepth until the end of
                                   our *owned* memory (before the mGhostDepth hyper slices start) */
                                  ptr + (mLayout.getSizesInMemory()[dimension]) * mLayout.stride(dimension),
                                  /* receive: in origin, including ghosts. */
                                  ptr);
      /* pointers: the same as above, but shifted by ghostDepth and ordering swapped. Yes. */
      mExchangeManager.exchangeDown(mGhostSubarrayMap.template getSubArray<T>(dimension), dimension,
                                    ptr + mGhostDepth * mLayout.stride(dimension),
                                    ptr + (mGhostDepth + mLayout.getSizesInMemory()[dimension]) *
                                              mLayout.stride(dimension));

      // BC-aware post-step: on boundary ranks, overwrite the wrap-around ghost slab with the BC
      // transform. Uses the same applyLocalBCAtDimDepth helper as Phase 2 / device path. Periodic
      // BC and non-boundary ranks are untouched.
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

  private:
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
