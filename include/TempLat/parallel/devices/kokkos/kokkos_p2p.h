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

#include <cstring>

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
