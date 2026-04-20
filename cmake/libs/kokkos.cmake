# ##############################################################################
# Get Kokkos
# ##############################################################################

message(
  "${BoldYellow}┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching Kokkos ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓${ColorReset}"
)

set(Kokkos_ENABLE_CUDA
    ${CUDA}
    CACHE BOOL "")
set(Kokkos_ENABLE_CUDA_CONSTEXPR
    ${CUDA}
    CACHE BOOL "")
set(Kokkos_ENABLE_HIP
    ${HIP}
    CACHE BOOL "")
set(Kokkos_ENABLE_OPENMP
    ${OPENMP}
    CACHE BOOL "")
set(Kokkos_ENABLE_THREADS
    ${PTHREADS}
    CACHE BOOL "")
set(Kokkos_ENABLE_SERIAL
    ${NOTHREADING}
    CACHE BOOL "")
set(Kokkos_ENABLE_TESTS
    OFF
    CACHE BOOL "")
set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION
    ON
    CACHE BOOL "")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(Kokkos_ENABLE_DEBUG
      ON
      CACHE BOOL "")
else()
  set(Kokkos_ENABLE_DEBUG
      OFF
      CACHE BOOL "")
endif()

FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  URL https://github.com/kokkos/kokkos/releases/download/5.1.0/kokkos-5.1.0.tar.gz
  URL_HASH
    SHA256=7bdbdfc88033ed7d940c7940ed8919e1f2b78a9656c69276beb76ad45c41ec4e
  ${SYSTEM_MARKER})
FetchContent_MakeAvailable(Kokkos)

message(
  "${BoldYellow}┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching Kokkos DONE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛${ColorReset}\n"
)
