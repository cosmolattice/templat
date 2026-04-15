# TempLat
*This repo is under active development; expect breaking changes.*
## Using TempLat in your project

*Minimal requirements:* 
- `CMake` version 3.16 or above
- `clang`, `g++` or another compiler with support for C++20

To include TempLat in your project, you can directly add it in your `CMakeLists.txt`, which will automatically download and have TempLat built as part of your project:

```cmake
cmake_minimum_required(VERSION 3.16)

project(my_project LANGUAGES CXX)

include(FetchContent)
set(FETCHCONTENT_QUIET ON)
FetchContent_Declare(
  TempLat
  GIT_REPOSITORY https://github.com/cosmolattice/templat.git
  GIT_TAG main)
FetchContent_MakeAvailable(TempLat)

# Add your executable
add_executable(my_executable main.cpp)
# Link your executable against TempLat
target_link_libraries(my_executable PRIVATE TempLat::TempLat)
```
As an example, your `main.cpp` can look like this:
```c++
#include "TempLat/session/sessionguard.h"

int main(int argc, char *argv[])
{
  TempLat::SessionGuard guard(argc, argv);

  sayMPI << "Hello, TempLat!\n";

  return 0;
};
```
And then you can build and run your project as usual with CMake:
```bash
mkdir build
cd build
cmake ..
make -j
./my_executable
```
Which should output
```
(0.000s) MPI Rank 0 - main.cpp:7 -->
Hello, TempLat!
```


### Choosing the device

By default, TempLat will attempt to detect what devices are available on a given machine.
In that case, it will first check for GPU support, checking first CUDA and then HIP. 
On the CPU side, it will first check for OpenMP support and then for C++ threads. If none of these are available, it will fall back to a serial implementation.
If you want to force the use of a specific device, you can do so by setting the appropriate flag when configuring the project with CMake.
The available options are `-DCUDA=ON`, `HIP=ON`, `-DOPENMP=ON`, `-DPTHREADS=ON` and `-DNOTHREADING=ON`. To force serial for example, you can do
```bash
cmake -DNOTHREADING=ON ..
```

#### CUDA

Using Nvidia GPUs requires that you have the CUDA toolkit, or the NVIDIA HPC SDK installed on your machine, or loaded as a module on your cluster. To compile for NVIDIA GPUs using CUDA, you can enable the CUDA backend by setting
```bash
cmake -DCUDA=ON ..
```
Specifying the architecture is optional for CUDA, as Kokkos can usually detect it correctly. However, if you want to specify it manually, you can do so by passing the appropriate flag to Kokkos as described in the section [Offline compilation (Kokkos)](#offline-compilation-kokkos) below.

#### HIP

Using AMD GPUs requires that you have the ROCm toolkit installed on your machine, or loaded as a module on your cluster.
To compile for AMD GPUs using HIP, you can enable the HIP backend by setting
```bash
export CXX=hipcc
cmake -DHIP=ON ..
```
Specifying the HIP compiler (`hipcc`) is necessary for the correct detection of the HIP toolkit.
Specifying the architecture is optional for HIP, as Kokkos can usually detect it correctly. However, if you want to specify it manually, you can do so by passing the appropriate flag to Kokkos as described in the section [Offline compilation (Kokkos)](#offline-compilation-kokkos) below.

#### cuFFTMp (multi-GPU FFTs via ParaFaFT)

When `PARAFAFT=ON` and `CUDA=ON`, ParaFaFT can additionally use NVIDIA's
[cuFFTMp](https://docs.nvidia.com/hpc-sdk/cufftmp/index.html) backend for
distributed multi-GPU FFTs. To enable it, pass `-DPARAFAFT_CUFFTMP=ON` at
configure time:
```bash
cmake -DMPI=ON -DCUDA=ON -DPARAFAFT=ON -DPARAFAFT_CUFFTMP=ON ..
```
The flag is forwarded unchanged to ParaFaFT. cuFFTMp and NVSHMEM must be
available on the system — they ship with the NVIDIA HPC SDK, or can be
located via the `CUFFTMP_HOME` / `NVSHMEM_HOME` environment variables (or
corresponding CMake cache variables). If ParaFaFT cannot find them the
configure step will fail with a clear error message.

### Offline compilation (Kokkos)

To compile an application to be run on a different architecture, you can directly pass the target architecture to Kokkos. For a list of supported architectures, see [the Kokkos documentation](https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html#gpu-architectures). For example, for an RTX 4070, you would pass 
```bash
cmake -DKokkos_ARCH_ADA89 ..
```
If no architecture is specified, Kokkos will attempt to detect the architecture of the machine. However, as compilation for GPU can take up to an hour, it is recommended on a cluster to queue the compilation separately on a node without a GPU, which requires offline compilation as described here.

### Device Providers

TempLat fully abstracts away the management and dispatching to actual devices, which is handled by device providers. The default device provider is Kokkos, which supports a wide range of devices and architectures. Currently, only Kokkos is implemented as a device provider.

### CMake Configuration Options

All custom CMake flags can be passed when configuring the user project, e.g. `cmake -DMPI=ON -DHDF5=ON ../`.

| Flag              | Description                        | Default                                |
| ----------------- | ---------------------------------- | -------------------------------------- |
| `PARAFAFT`        | ParaFaft support for parallel FFTs | `OFF`                                  |
| `MPI`             | MPI support                        | `OFF`                                  |
| `HDF5`            | HDF5 support                       | `OFF`                                  |
| `TEMPLAT_TEST`   | Enable TempLat's tests             | `OFF`                                  |
| `DEVICE_PROVIDER` | Backend for parallelization        | `Kokkos`                               |
| `CUDA`            | CUDA support for NVIDIA GPUs       | `OFF`                                  |
| `HIP`             | HIP support for AMD GPUs           | `OFF`                                  |
| `OPENMP`          | OpenMP CPU parallelization         | `OFF`                                  |
| `PTHREADS`        | C++ threads CPU parallelization    | `OFF`                                  |
| `NOTHREADING`     | No parallelization                 | `OFF`                                  |
| `NATIVE`          | Pass `--march=native` to compiler  | `ON` (non-macOS), `OFF` (macOS)        |
| `KOKKOSFFT`       | KokkosFFT for single-node GPU FFTs | `ON` when CUDA/HIP enabled, else `OFF` |
| `PARAFAFT_CUFFTMP`| cuFFTMp backend inside ParaFaFT    | `OFF`                                  |

### Runtime Environment Variables

| Variable          | Description                                                                                                                                                                                                                                                                      |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `GPU_NOCONSTRAIN` | Set to `1` to allow more MPI ranks per node than GPUs (oversubscription, e.g. via CUDA MPS). By default the session aborts if `local_ranks > num_devices`. Device assignment wraps via `shmrank % num_devices`; same-GPU peers use CUDA IPC directly with no MPI on the halo path. |
