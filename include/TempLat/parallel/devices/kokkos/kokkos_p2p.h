#ifndef TEMPLAT_PARALLEL_KOKKOS_P2P_H
#define TEMPLAT_PARALLEL_KOKKOS_P2P_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2026

#include "TempLat/util/exception.h"

#if (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)) && defined(HAVE_MPI)

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_runtime.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

namespace TempLat::device_kokkos::p2p
{
  MakeException(GpuP2PException);

  // ============================================================
  // Backend dispatch: thin wrappers mapping CUDA <-> HIP uniformly
  // ============================================================

#if defined(KOKKOS_ENABLE_CUDA)

  using IpcMemHandle_t = cudaIpcMemHandle_t;
  static constexpr size_t IpcHandleSize = sizeof(cudaIpcMemHandle_t);

  inline int getDeviceId()
  {
    int dev;
    auto err = cudaGetDevice(&dev);
    if (err != cudaSuccess)
      throw GpuP2PException("cudaGetDevice failed: ", cudaGetErrorString(err));
    return dev;
  }

  inline bool canAccessPeer(int srcDevice, int dstDevice)
  {
    int canAccess = 0;
    cudaDeviceCanAccessPeer(&canAccess, srcDevice, dstDevice);
    return canAccess != 0;
  }

  inline void enablePeerAccess(int peerDevice)
  {
    cudaError_t err = cudaDeviceEnablePeerAccess(peerDevice, 0);
    if (err != cudaSuccess && err != cudaErrorPeerAccessAlreadyEnabled)
      throw GpuP2PException("cudaDeviceEnablePeerAccess failed for device ", peerDevice, ": ",
                            cudaGetErrorString(err));
  }

  inline void ipcGetHandle(void *devPtr, void *handle)
  {
    auto err = cudaIpcGetMemHandle(reinterpret_cast<cudaIpcMemHandle_t *>(handle), devPtr);
    if (err != cudaSuccess)
      throw GpuP2PException("cudaIpcGetMemHandle failed: ", cudaGetErrorString(err));
  }

  inline void *ipcOpenHandle(const void *handle)
  {
    void *ptr;
    auto err = cudaIpcOpenMemHandle(&ptr, *reinterpret_cast<const cudaIpcMemHandle_t *>(handle),
                                    cudaIpcMemLazyEnablePeerAccess);
    if (err != cudaSuccess)
      throw GpuP2PException("cudaIpcOpenMemHandle failed: ", cudaGetErrorString(err));
    return ptr;
  }

  inline void ipcCloseHandle(void *ptr) { cudaIpcCloseMemHandle(ptr); }

  inline void memcpyAsync(void *dst, const void *src, size_t bytes)
  {
    auto err = cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDefault, 0);
    if (err != cudaSuccess)
      throw GpuP2PException("cudaMemcpyAsync failed: ", cudaGetErrorString(err));
  }

  inline void streamSynchronize()
  {
    auto err = cudaStreamSynchronize(0);
    if (err != cudaSuccess)
      throw GpuP2PException("cudaStreamSynchronize failed: ", cudaGetErrorString(err));
  }

#elif defined(KOKKOS_ENABLE_HIP)

  using IpcMemHandle_t = hipIpcMemHandle_t;
  static constexpr size_t IpcHandleSize = sizeof(hipIpcMemHandle_t);

  inline int getDeviceId()
  {
    int dev;
    auto err = hipGetDevice(&dev);
    if (err != hipSuccess)
      throw GpuP2PException("hipGetDevice failed: ", hipGetErrorString(err));
    return dev;
  }

  inline bool canAccessPeer(int srcDevice, int dstDevice)
  {
    int canAccess = 0;
    hipDeviceCanAccessPeer(&canAccess, srcDevice, dstDevice);
    return canAccess != 0;
  }

  inline void enablePeerAccess(int peerDevice)
  {
    hipError_t err = hipDeviceEnablePeerAccess(peerDevice, 0);
    if (err != hipSuccess && err != hipErrorPeerAccessAlreadyEnabled)
      throw GpuP2PException("hipDeviceEnablePeerAccess failed for device ", peerDevice, ": ",
                            hipGetErrorString(err));
  }

  inline void ipcGetHandle(void *devPtr, void *handle)
  {
    auto err = hipIpcGetMemHandle(reinterpret_cast<hipIpcMemHandle_t *>(handle), devPtr);
    if (err != hipSuccess)
      throw GpuP2PException("hipIpcGetMemHandle failed: ", hipGetErrorString(err));
  }

  inline void *ipcOpenHandle(const void *handle)
  {
    void *ptr;
    auto err =
        hipIpcOpenMemHandle(&ptr, *reinterpret_cast<const hipIpcMemHandle_t *>(handle), hipIpcMemLazyEnablePeerAccess);
    if (err != hipSuccess)
      throw GpuP2PException("hipIpcOpenMemHandle failed: ", hipGetErrorString(err));
    return ptr;
  }

  inline void ipcCloseHandle(void *ptr) { hipIpcCloseMemHandle(ptr); }

  inline void memcpyAsync(void *dst, const void *src, size_t bytes)
  {
    auto err = hipMemcpyAsync(dst, src, bytes, hipMemcpyDefault, 0);
    if (err != hipSuccess)
      throw GpuP2PException("hipMemcpyAsync failed: ", hipGetErrorString(err));
  }

  inline void streamSynchronize()
  {
    auto err = hipStreamSynchronize(0);
    if (err != hipSuccess)
      throw GpuP2PException("hipStreamSynchronize failed: ", hipGetErrorString(err));
  }

#endif

  // ============================================================
  // Raw GPU memory allocation (bypasses Kokkos header)
  // ============================================================

  // Kokkos::View.data() is offset from the cudaMalloc base by a SharedAllocationHeader
  // (128 bytes). cudaIpcGetMemHandle/cudaIpcOpenMemHandle operate on the allocation base.
  // To avoid this offset issue, IPC-exported buffers are allocated with raw cudaMalloc/hipMalloc
  // and wrapped as unmanaged Kokkos views. This gives us clean base pointers for IPC.

  inline void *rawDeviceMalloc(size_t bytes)
  {
    void *ptr = nullptr;
#if defined(KOKKOS_ENABLE_CUDA)
    auto err = cudaMalloc(&ptr, bytes);
    if (err != cudaSuccess)
      throw GpuP2PException("cudaMalloc failed for IPC buffer (", bytes, " bytes): ", cudaGetErrorString(err));
#elif defined(KOKKOS_ENABLE_HIP)
    auto err = hipMalloc(&ptr, bytes);
    if (err != hipSuccess)
      throw GpuP2PException("hipMalloc failed for IPC buffer (", bytes, " bytes): ", hipGetErrorString(err));
#endif
    return ptr;
  }

  inline void rawDeviceFree(void *ptr)
  {
    if (ptr == nullptr) return;
#if defined(KOKKOS_ENABLE_CUDA)
    cudaFree(ptr);
#elif defined(KOKKOS_ENABLE_HIP)
    hipFree(ptr);
#endif
  }

  // ============================================================
  // Link type detection: NVLink/xGMI (full-duplex) vs PCIe
  // ============================================================

  inline bool isFullDuplexLink(int deviceA, int deviceB)
  {
#if defined(KOKKOS_ENABLE_CUDA)
    // cudaDevP2PAttrCudaArrayAccessSupported is only true over NVLink/NVSwitch, not PCIe.
    int value = 0;
    cudaDeviceGetP2PAttribute(&value, cudaDevP2PAttrCudaArrayAccessSupported, deviceA, deviceB);
    return value != 0;
#elif defined(KOKKOS_ENABLE_HIP)
    uint32_t linktype = 0, hopcount = 0;
    auto err = hipExtGetLinkTypeAndHopCount(deviceA, deviceB, &linktype, &hopcount);
    if (err != hipSuccess) return false;
    // HSA_AMD_LINK_INFO_TYPE_XGMI = 4 (AMD Infinity Fabric, full-duplex)
    return linktype == 4;
#endif
  }

  // ============================================================
  // IPC handle packet — POD struct sent via MPI_BYTE
  // ============================================================

  struct IpcHandlePacket {
    char handle[IpcHandleSize];
    int deviceId;
    uint64_t version;
  };

} // namespace TempLat::device_kokkos::p2p

#endif // (KOKKOS_ENABLE_CUDA || KOKKOS_ENABLE_HIP) && HAVE_MPI

#endif
