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

function(templat_define_options)
  option(CUDA "Enable CUDA support" OFF)
  option(HIP "Enable HIP support" OFF)
  option(OpenMP "Enable OpenMP support" OFF)
  option(Threads "Enable Threads support" OFF)
  option(Serial "Enable Serial support" OFF)

  set(DEVICE_PROVIDER
      "Kokkos"
      CACHE STRING "Device provider to use (e.g. Kokkos)")
  set(KOKKOSFFT
      ON
      CACHE STRING "Enable KokkosFFT")

  option(HDF5 "Use HDF5 for I/O" OFF)
  option(MPI "Use MPI for parallelization" OFF)
  option(PARAFAFT "Use ParaFaft for parallel FFTs" OFF)
  option(TEMPLAT_TESTS "Build TempLat tests" OFF)
endfunction()

templat_define_options()

# ##############################################################################
# Create the TempLat (interface) library
# ##############################################################################

add_library(TempLat INTERFACE)
add_library(TempLat::TempLat ALIAS TempLat)

# ##############################################################################
# Get all dependencies
# ##############################################################################

set(THREADS_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)
target_link_libraries(TempLat INTERFACE Threads::Threads)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/fftw.cmake)
target_link_libraries(TempLat INTERFACE ${FFTW_LIBRARIES})

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/device.cmake)
target_link_device(TempLat)

if(MPI)
  target_link_libraries(TempLat INTERFACE ${MPI_CXX_LIBRARIES})
endif()

if(HDF5)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5.cmake)
  target_link_libraries(TempLat INTERFACE ${HDF5_LIBRARIES})
endif()

if(MPI AND PARAFAFT)
  target_link_libraries(TempLat INTERFACE parafaft::parafaft)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/parafaf.cmake)
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
