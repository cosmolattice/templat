# ##############################################################################
# Find out which backends are available
#
# If nothing is selected, we choose CUDA > HIP > OpenMP > Threads > Serial
# ##############################################################################

include(CheckLanguage)
message("")

if(NOT CUDA
   AND NOT HIP
   AND NOT OpenMP
   AND NOT Threads
   AND NOT Serial)
  message(
    STATUS
      "${BoldCyan}━━━━━━━━━━━━━━━━━━━ No device specified, trying to auto-detect ━━━━━━━━━━━━━━━━━━━${ColorReset}"
  )
  set(CUDA ON)
  set(HIP ON)
  set(OpenMP ON)
  set(Threads ON)
  set(Serial ON)
endif()

if(CUDA)
  # Let's see if we have a CUDA compiler
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    set(CUDA ON)
    set(HIP OFF)
    set(OpenMP ON)
    set(Threads OFF)
    set(Serial OFF)
  else()
    set(CUDA OFF)
  endif()
endif()

if(NOT CUDA AND HIP)
  # There is a problem with Kokkos + HIP, where HIP's cmake sets offload-arch
  # and Kokkos does the same, which leads to a conflict. We purge the info set
  # by HIP's cmake here.
  if(HIP)
    set(GPU_BUILD_TARGETS
        ""
        CACHE STRING "" FORCE)
    set(GPU_TARGETS
        ""
        CACHE STRING "" FORCE)
  endif()

  # Let's see if we have a HIP compiler
  check_language(HIP)
  if(CMAKE_HIP_COMPILER)
    set(CUDA OFF)
    set(HIP ON)
    set(OpenMP ON)
    set(Threads OFF)
    set(Serial OFF)
    # For HIP, we also need to specify the platform (amd or nvidia). Kokkos only
    # supports amd as an option here, though.
    set(HIP_PLATFORM
        "amd"
        CACHE STRING "Set the HIP platform (amd or nvidia)")

  else()
    set(HIP OFF)
  endif()
endif()

if(NOT CUDA
   AND NOT HIP
   AND OpenMP)
  check_language(OpenMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    set(OpenMP ON)
    set(Threads OFF)
    set(Serial OFF)
  else()
    set(OpenMP OFF)
  endif()
endif()

if(NOT CUDA
   AND NOT HIP
   AND NOT OpenMP
   AND Serial)
  set(Serial ON)
  set(CUDA OFF)
  set(HIP OFF)
  set(OpenMP OFF)
  set(Threads OFF)
endif()

if(NOT CUDA
   AND NOT HIP
   AND NOT OpenMP
   AND NOT Serial
   AND Threads)
  set(CUDA OFF)
  set(HIP OFF)
  set(OpenMP OFF)
  set(Serial OFF)
  set(Threads ON)
endif()

if(NOT CUDA
   AND NOT HIP
   AND NOT OpenMP
   AND NOT Threads
   AND NOT Serial)
  message(
    FATAL_ERROR
      "No valid device configuration found. Please specify at least one of CUDA, HIP, OpenMP, Threads or Serial."
  )
endif()

message(
  STATUS
    "Device configuration: \n    CUDA: ${CUDA} \n    HIP: ${HIP} \n    OpenMP: ${OpenMP} \n    Threads: ${Threads} \n    Serial: ${Serial}"
)
