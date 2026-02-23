cmake_minimum_required(VERSION 3.16)
project(
  TempLat
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/common.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/title.cmake)

# Options and configuration

option(HDF5 "Use HDF5 for I/O" OFF)
option(PARAFAFT "Use ParaFaft for parallel FFTs" OFF)
set(DEVICE_PROVIDER
    "Kokkos"
    CACHE STRING "Device provider to use (e.g. Kokkos)")

# ---------- Get required dependencies ----------
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/fftw.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/device.cmake)

if(HDF5)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5.cmake)
endif()

if(PARAFAFT)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/parafaf.cmake)
endif()

# ---------- Header-only interface library ----------
add_library(TempLat INTERFACE)
add_library(TempLat::TempLat ALIAS TempLat)

target_include_directories(
  TempLat
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/include>
            $<INSTALL_INTERFACE:include>)

target_link_libraries(TempLat INTERFACE Threads::Threads)
target_link_libraries(TempLat INTERFACE ${FFTW_LIBRARIES})
target_link_device(TempLat)
if(MPI)
  target_link_libraries(TempLat INTERFACE ${MPI_CXX_LIBRARIES})
endif()
if(PARAFAFT AND MPI)
  target_link_libraries(TempLat INTERFACE parafaft::parafaft)
endif()
if(HDF5)
  target_link_libraries(TempLat INTERFACE ${HDF5_LIBRARIES})
endif()

# ---------- Tests (off by default, e.g. when consumed via FetchContent)
# ----------
option(TEMPLAT_BUILD_TESTS "Build TempLat tests" OFF)

if(TEMPLAT_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
