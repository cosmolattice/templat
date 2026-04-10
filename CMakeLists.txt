# ##############################################################################
# CMake configuration for TempLat, a modern C++ template library for lattice
# computations.
#
# CMake Options: - PARAFAFT: compile TempLat with ParaFaft support for parallel
# FFTs. Default: OFF - MPI: compile TempLat with MPI support. Default: OFF -
# HDF5: compile TempLat with HDF5 support. Default: OFF - TEMPLAT_TESTS: enable
# TempLat's tests. Default: OFF - DEVICE_PROVIDER: switches the backend for
# parallelization. Default: Kokkos - CUDA: compile with CUDA support for GPU
# parallelization. Default: OFF - HIP: compile with HIP support for GPU
# parallelization. Default: OFF - OPENMP: compile with OpenMP support for CPU
# parallelization. Default: OFF - PTHREADS: compile with C++ threads support for
# CPU parallelization. Default: OFF - NOTHREADING: compile without
# parallelization. Default: OFF By default, TempLat will try to auto-detect the
# system's capabilities and choose the backend. - FLOAT: compile with
# single-precision FFTW support. Default: OFF - NATIVE: pass --march=native to
# the compiler for better performance. Default: ON on non-Darwin systems, OFF on
# Darwin (macOS) - KOKKOSFFT: compile TempLat with KokkosFFT support for
# single-node GPU FFTs. Default: ON when CUDA or HIP support is enabled, OFF
# otherwise
#
# ##############################################################################

cmake_minimum_required(VERSION 3.16)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/common.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/title.cmake)

if((NOT DEFINED CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL ""))
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Build type" FORCE)
  message(STATUS "No build type specified, defaulting to Release.")
endif()

project(
  TempLat
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ##############################################################################
# Options that can be set by the user
# ##############################################################################

option(CUDA "Enable CUDA support" OFF)
option(HIP "Enable HIP support" OFF)
option(OPENMP "Enable OpenMP support" OFF)
option(PTHREADS "Enable pthreads support" OFF)
option(NOTHREADING "Disable threading" OFF)

option(FLOAT "Enable single-precision (float) FFTW support" OFF)

set(DEVICE_PROVIDER
    "Kokkos"
    CACHE STRING "Device provider to use (e.g. Kokkos)")
# KokkosFFT is set below based on the device backend!

option(AUTOBUILD_FFTW
       "Automatically build FFTW from source if not found on the system" OFF)
option(AUTOBUILD_HDF5
       "Automatically build HDF5 from source if not found on the system" OFF)

option(HDF5 "Use HDF5 for I/O" OFF)
option(MPI "Use MPI for parallelization" OFF)
option(PARAFAFT "Use ParaFaft for parallel FFTs" OFF)
option(TEMPLAT_TESTS "Build TempLat tests" OFF)
option(TEMPLAT_BENCH "Build TempLat benchmarks" OFF)

# ##############################################################################
# Compatibility with CL1.0
# ##############################################################################

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compatibility.cmake)

# ##############################################################################
# Create the TempLat (interface) library
# ##############################################################################

add_library(TempLat INTERFACE)
add_library(TempLat::TempLat ALIAS TempLat)

# ##############################################################################
# Get all dependencies
# ##############################################################################

# Start with understanding which backends are available and should be used
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/device.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/device_provider.cmake)

# Now that we know which device provider is being used, we can set KOKKOSFFT
# accordingly
if(CUDA OR HIP)
  set(KOKKOSFFT
      ON
      CACHE STRING "Enable KokkosFFT")
else()
  set(KOKKOSFFT
      OFF
      CACHE STRING "Enable KokkosFFT")
endif()

# We use threads for FFTW, as well as for detecting hardware concurrency
set(THREADS_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)
target_link_libraries(TempLat INTERFACE Threads::Threads)

# We need to force on float support for FFTW if KokkosFFT is enabled, since
# KokkosFFT requires it.
if(KOKKOSFFT)
  if(NOT FLOAT)
    set(FLOAT
        ON
        CACHE
          BOOL
          "Enable single-precision (float) FFTW support (required by KokkosFFT)"
          FORCE)
  endif()
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/fftw.cmake)
target_link_libraries(TempLat INTERFACE ${FFTW_LIBRARIES})
target_link_directories(TempLat INTERFACE ${FFTW_LIB_DIR})
target_include_directories(TempLat SYSTEM BEFORE INTERFACE ${FFTW_INCLUDES})

if(MPI)
  find_package(
    MPI
    COMPONENTS CXX
    REQUIRED)
  target_include_directories(TempLat SYSTEM BEFORE
                             INTERFACE ${MPI_CXX_INCLUDE_DIRS})
  target_link_libraries(TempLat INTERFACE ${MPI_CXX_LIBRARIES})
  target_compile_definitions(TempLat INTERFACE HAVE_MPI)
  message(
    STATUS "MPI enabled, including MPI headers from ${MPI_CXX_INCLUDE_DIRS}")

  # Force ParaFaFT on if MPI is enabled, since ParaFaFT is just better than FFTW
  # and there's no reason not to use it if MPI is available.
  if(NOT PARAFAFT)
    set(PARAFAFT
        ON
        CACHE BOOL "Enable ParaFaFT support (recommended if MPI is enabled)"
              FORCE)
  endif()
endif()

if(KOKKOSFFT)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/kokkos-fft.cmake)
  target_link_libraries(TempLat INTERFACE KokkosFFT::fft)
  target_compile_definitions(TempLat INTERFACE HAVE_KOKKOSFFT)
endif()

if(HDF5)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5.cmake)
  target_link_libraries(TempLat INTERFACE ${HDF5_LIBRARIES})
  target_compile_definitions(TempLat INTERFACE HAVE_HDF5)
endif()

if(MPI AND PARAFAFT)
  target_link_libraries(TempLat INTERFACE parafaft::parafaft)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/parafaft.cmake)
  target_compile_definitions(TempLat INTERFACE HAVE_PARAFAFT)
elseif(PARAFAFT)
  message(
    WARNING
      "ParaFaFT requires MPI. ParaFaFT will not be enabled since MPI was not enabled."
  )
endif()

target_include_directories(
  TempLat INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_include_directories(
  TempLat SYSTEM
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/include>
            $<INSTALL_INTERFACE:include>)

# ##############################################################################
# System-specific compiler options.
# ##############################################################################

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/system.cmake)

# ##############################################################################
# Tests (Off by default, see above)
# ##############################################################################

if(TEMPLAT_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(TEMPLAT_BENCH)
  add_subdirectory(benchmarks)
endif()
