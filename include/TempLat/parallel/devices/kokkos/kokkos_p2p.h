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
#include <cstdio>
#include <dlfcn.h>
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

#if defined(KOKKOS_ENABLE_CUDA)
  // CUDA runtime P2P attributes (cudaDevP2PAttrCudaArrayAccessSupported, etc.) are
  // unreliable for distinguishing NVLink from PCIe. Query NVML's NVLink port state
  // via dlopen — libnvidia-ml.so.1 is always present with the NVIDIA driver.
  inline bool isFullDuplexLink(int cudaDeviceA, int cudaDeviceB)
  {
    // Get PCI addresses from CUDA to match against NVML devices
    cudaDeviceProp propA{}, propB{};
    if (cudaGetDeviceProperties(&propA, cudaDeviceA) != cudaSuccess) return false;
    if (cudaGetDeviceProperties(&propB, cudaDeviceB) != cudaSuccess) return false;

    void *lib = dlopen("libnvidia-ml.so.1", RTLD_LAZY);
    if (!lib) return false; // NVML unavailable → assume PCIe (safe default)

    // NVML function signatures (avoid header dependency)
    using InitFn = unsigned int (*)();
    using ShutdownFn = unsigned int (*)();
    using HandleByPciFn = unsigned int (*)(const char *, void **);
    using NvLinkStateFn = unsigned int (*)(void *, unsigned int, unsigned int *);
    using NvLinkRemotePciFn = unsigned int (*)(void *, unsigned int, void *);
    using GetPciInfoFn = unsigned int (*)(void *, void *);

    auto fnInit = reinterpret_cast<InitFn>(dlsym(lib, "nvmlInit_v2"));
    auto fnShutdown = reinterpret_cast<ShutdownFn>(dlsym(lib, "nvmlShutdown"));
    auto fnHandleByPci = reinterpret_cast<HandleByPciFn>(dlsym(lib, "nvmlDeviceGetHandleByPciBusId_v2"));
    auto fnNvLinkState = reinterpret_cast<NvLinkStateFn>(dlsym(lib, "nvmlDeviceGetNvLinkState"));
    auto fnNvLinkRemotePci = reinterpret_cast<NvLinkRemotePciFn>(dlsym(lib, "nvmlDeviceGetNvLinkRemotePciInfo_v2"));
    auto fnGetPciInfo = reinterpret_cast<GetPciInfoFn>(dlsym(lib, "nvmlDeviceGetPciInfo_v3"));

    if (!fnInit || !fnShutdown || !fnHandleByPci || !fnNvLinkState || !fnNvLinkRemotePci || !fnGetPciInfo) {
      dlclose(lib);
      return false;
    }

    bool found = false;
    if (fnInit() == 0) { // NVML_SUCCESS = 0
      // Look up NVML device handles by PCI bus ID (handles CUDA_VISIBLE_DEVICES reordering)
      char busIdA[32], busIdB[32];
      snprintf(busIdA, sizeof(busIdA), "%08x:%02x:%02x.0", propA.pciDomainID, propA.pciBusID, propA.pciDeviceID);
      snprintf(busIdB, sizeof(busIdB), "%08x:%02x:%02x.0", propB.pciDomainID, propB.pciBusID, propB.pciDeviceID);

      void *devA = nullptr, *devB = nullptr;
      if (fnHandleByPci(busIdA, &devA) == 0 && fnHandleByPci(busIdB, &devB) == 0) {
        // Get device B's PCI info for comparison with NVLink remote endpoints
        // nvmlPciInfo_t layout: char busIdLegacy[16], uint domain, uint bus, uint device, ...
        struct { char legacy[16]; unsigned int domain, bus, device; } pciB{};
        fnGetPciInfo(devB, &pciB);

        // Enumerate NVLink ports on device A (up to 18 on recent hardware)
        for (unsigned int link = 0; link < 18 && !found; ++link) {
          unsigned int state = 0; // nvmlEnableState_t
          if (fnNvLinkState(devA, link, &state) != 0 || state != 1) // NVML_FEATURE_ENABLED = 1
            continue;
          struct { char legacy[16]; unsigned int domain, bus, device; } remotePci{};
          if (fnNvLinkRemotePci(devA, link, &remotePci) != 0)
            continue;
          if (remotePci.domain == pciB.domain && remotePci.bus == pciB.bus && remotePci.device == pciB.device)
            found = true;
        }
      }
      fnShutdown();
    }
    dlclose(lib);
    return found;
  }

#elif defined(KOKKOS_ENABLE_HIP)
  inline bool isFullDuplexLink(int deviceA, int deviceB)
  {
    uint32_t linktype = 0, hopcount = 0;
    auto err = hipExtGetLinkTypeAndHopCount(deviceA, deviceB, &linktype, &hopcount);
    if (err != hipSuccess) return false;
    // HSA_AMD_LINK_INFO_TYPE_XGMI = 4 (AMD Infinity Fabric, full-duplex)
    return linktype == 4;
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
