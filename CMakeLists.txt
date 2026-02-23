cmake_minimum_required(VERSION 3.16)
project(TempLat VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------- Header-only interface library ----------
add_library(TempLat INTERFACE)
add_library(TempLat::TempLat ALIAS TempLat)

target_include_directories(TempLat INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/include>
    $<INSTALL_INTERFACE:include>
)

# ---------- Tests (off by default, e.g. when consumed via FetchContent) ----------
option(TEMPLAT_BUILD_TESTS "Build TempLat tests" OFF)

if(TEMPLAT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
