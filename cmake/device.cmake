# Device detection + Kokkos fetch for standalone TempLat test builds.
# NOT loaded when TempLat is consumed via FetchContent (TEMPLAT_BUILD_TESTS=OFF).
include(CheckLanguage)

option(CUDA "Enable CUDA support" OFF)
option(HIP "Enable HIP support" OFF)
option(OpenMP "Enable OpenMP support" OFF)
option(Threads "Enable Threads support" OFF)
option(Serial "Enable Serial support" OFF)

# Auto-detect if nothing specified
if(NOT CUDA AND NOT HIP AND NOT OpenMP AND NOT Threads AND NOT Serial)
    set(CUDA ON)
    set(HIP ON)
    set(OpenMP ON)
    set(Threads ON)
    set(Serial ON)
endif()

if(CUDA)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        set(HIP OFF)
        set(OpenMP ON)
        set(Threads OFF)
        set(Serial OFF)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEVICE_CUDA")
    else()
        set(CUDA OFF)
    endif()
endif()

if(NOT CUDA AND HIP)
    # Purge HIP's cmake arch settings to avoid conflict with Kokkos
    if(HIP)
        set(GPU_BUILD_TARGETS "" CACHE STRING "" FORCE)
        set(GPU_TARGETS "" CACHE STRING "" FORCE)
    endif()

    check_language(HIP)
    if(CMAKE_HIP_COMPILER)
        set(OpenMP ON)
        set(Threads OFF)
        set(Serial OFF)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEVICE_HIP")
        set(HIP_PLATFORM "amd" CACHE STRING "Set the HIP platform (amd or nvidia)")
    else()
        set(HIP OFF)
    endif()
endif()

if(NOT CUDA AND NOT HIP AND OpenMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        set(Threads OFF)
        set(Serial OFF)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEVICE_CPU")
    else()
        set(OpenMP OFF)
    endif()
endif()

if(NOT CUDA AND NOT HIP AND NOT OpenMP AND Threads)
    set(Serial OFF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEVICE_CPU")
endif()

if(NOT CUDA AND NOT HIP AND NOT OpenMP AND NOT Threads AND Serial)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEVICE_CPU")
endif()

if(NOT CUDA AND NOT HIP AND NOT OpenMP AND NOT Threads AND NOT Serial)
    message(FATAL_ERROR "No valid device configuration found.")
endif()

message(STATUS "TempLat device: CUDA=${CUDA} HIP=${HIP} OpenMP=${OpenMP} Threads=${Threads} Serial=${Serial}")

# ---------- Fetch Kokkos ----------
if(NOT TARGET Kokkos::kokkos)
    message(STATUS "Fetching Kokkos for TempLat tests...")
    set(Kokkos_ENABLE_CUDA ${CUDA} CACHE BOOL "")
    set(Kokkos_ENABLE_CUDA_CONSTEXPR ${CUDA} CACHE BOOL "")
    set(Kokkos_ENABLE_HIP ${HIP} CACHE BOOL "")
    set(Kokkos_ENABLE_OPENMP ${OpenMP} CACHE BOOL "")
    set(Kokkos_ENABLE_THREADS ${Threads} CACHE BOOL "")
    set(Kokkos_ENABLE_SERIAL ${Serial} CACHE BOOL "")
    set(Kokkos_ENABLE_TESTS OFF CACHE BOOL "")
    set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION ON CACHE BOOL "")

    include(FetchContent)
    FetchContent_Declare(Kokkos
        DOWNLOAD_EXTRACT_TIMESTAMP FALSE
        URL https://github.com/kokkos/kokkos/releases/download/5.0.2/kokkos-5.0.2.tar.gz
        URL_HASH SHA256=188817bb452ca805ee8701f1c5adbbb4fb83dc8d1c50624566a18a719ba0fa5e
        SYSTEM)
    FetchContent_MakeAvailable(Kokkos)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEVICE_KOKKOS")

# ---------- Optional: KokkosFFT ----------
option(KOKKOSFFT "Enable KokkosFFT support" ON)
if(KOKKOSFFT AND NOT TARGET KokkosFFT::fft)
    message(STATUS "Fetching KokkosFFT for TempLat tests...")
    include(FetchContent)
    FetchContent_Declare(KokkosFFT
        DOWNLOAD_EXTRACT_TIMESTAMP FALSE
        URL https://github.com/kokkos/kokkos-fft/archive/refs/tags/v1.0.0.zip
        URL_HASH SHA256=80e9c1abdf71df2342ae713c845ba2aeabf1c1a0c2e116795534a8e67734c177
        SYSTEM)
    FetchContent_MakeAvailable(KokkosFFT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_KOKKOSFFT")
    set(KOKKOSFFT ON)
endif()

# ---------- Helper function to link device libraries ----------
function(target_link_device target)
    target_link_libraries(${target} PUBLIC Kokkos::kokkos)
    if(KOKKOSFFT AND TARGET KokkosFFT::fft)
        target_link_libraries(${target} PUBLIC KokkosFFT::fft)
    endif()
endfunction()
