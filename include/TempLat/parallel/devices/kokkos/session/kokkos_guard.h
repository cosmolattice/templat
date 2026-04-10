#ifndef TEMPLAT_PARALLEL_KOKKOS_SESSION_KOKKOSGUARD_H
#define TEMPLAT_PARALLEL_KOKKOS_SESSION_KOKKOSGUARD_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include "TempLat/util/exception.h"
#include "TempLat/parallel/devices/kokkos/kokkos.h"
#include "TempLat/parallel/threadsettings.h"
#include "TempLat/util/log/saycomplete.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace TempLat::device_kokkos
{
  MakeException(KokkosDeviceGuardInstantiationException);

  /**
   * @brief A class to manage and initialize the Kokkos runtime environment and its dependencies.
   *
   * The DeviceGuard class ensures proper initialization and finalization of the Kokkos runtime
   * environment and associated libraries such as MPI. It is designed to be instantiated exactly
   * once per process, typically at the beginning of the main function. The class handles thread
   * safety, MPI threading level verification, and manages a local sub-communicator for node-local
   * operations.
   *
   * Proper use of this class prevents multiple unwanted instances of Kokkos initialization and
   * ensures that all resources are cleaned up upon program termination.
   */
  class DeviceGuard
  {
  public:
    // Put public methods here. These should change very little over time.

    DeviceGuard(int argc, char **argv, bool verbose = false)
        : instanceProtectionKey(InstanceCounter(1)), mVerbose(verbose)
    {
      auto threadSettings = ThreadSettings::getInstance();

      Kokkos::InitializationSettings kokkos_settings;
      kokkos_settings.set_print_configuration(this->mVerbose);
      kokkos_settings.set_num_threads(threadSettings.getMaxThreadCount());

// We need to do load-balancing here, if we are using GPU + MPI
#ifdef HAVE_MPI
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      int num_devices = Kokkos::num_devices();

      // First, make an MPI group for the local machine.
      MPI_Comm shmcomm;
      MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
      sShmComm = shmcomm;
      int shmrank, shmsize;
      MPI_Comm_rank(shmcomm, &shmrank);
      MPI_Comm_size(shmcomm, &shmsize);

      int gpuNoConstrain = 0;
      try {
        if (const char *env_p = std::getenv("GPU_NOCONSTRAIN")) gpuNoConstrain = std::stoi(env_p);
      } catch (...) {
        throw KokkosDeviceGuardInstantiationException(
            "Error parsing GPU_NOCONSTRAIN environment variable. Expected an integer. Got: ",
            std::getenv("GPU_NOCONSTRAIN"));
      }

      // Check if there are more processes on this node than devices.
      if (num_devices < shmsize && gpuNoConstrain == 0)
        throw KokkosDeviceGuardInstantiationException(
            "There are more MPI processes on this node than available GPU devices. This will lead to "
            "problems. Number of devices: ",
            num_devices, ", number of processes on this node: ", shmsize);

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      // Assign devices to processes in a round-robin fashion.
      sDeviceId = shmrank % num_devices;
      kokkos_settings.set_device_id(sDeviceId);
      std::cout << "Global process rank " << rank << " (local rank " << shmrank << " on this node) assigned to device "
                << sDeviceId << std::endl;
#endif
#endif

      Kokkos::initialize(kokkos_settings);
    }

    ~DeviceGuard() { Kokkos::finalize(); }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    int instanceProtectionKey;
    bool mVerbose;
    /** @brief A sub-group of all the processes, which is local on the same node. This can be used to compute the best
     number of threads per process. By default we are greedy: use all the cpu power on a node. Set arguments upon
     launch to avoid that. (Not implemented yet.) */

    static inline int InstanceCounter(int delta = 0)
    {
      static int counter = 0;
      counter += delta;
      if (counter > 1)
        throw KokkosDeviceGuardInstantiationException(
            "Per process, Kokkos' DeviceGuard can be instantiated only once. This "
            "should be done in `int main()`. This is wrong. Instances:",
            counter);
      return counter;
    }

  public:
    static int GetInstanceCount() { return InstanceCounter(); }

#ifdef HAVE_MPI
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    static MPI_Comm getShmComm() { return sShmComm; }
    static int getDeviceId() { return sDeviceId; }

  private:
    static inline MPI_Comm sShmComm = MPI_COMM_NULL;
    static inline int sDeviceId = -1;
#else
  public:
    static MPI_Comm getShmComm() { return MPI_COMM_NULL; }
    static int getDeviceId() { return 0; }
#endif
#else
  public:
    static int getDeviceId() { return 0; }
#endif
  };
} // namespace TempLat::device_kokkos

#endif
