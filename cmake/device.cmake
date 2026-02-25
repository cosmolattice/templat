# ##############################################################################
# Find out which backends are available
#
# Backends are split into two independent categories:
#   GPU: CUDA, HIP               (at most one;  priority: CUDA > HIP)
#   CPU: OpenMP, Threads, Serial  (exactly one;  priority: OpenMP > Threads > Serial)
#
# If nothing is specified, both categories are auto-detected.
# If the user specifies only a GPU backend, the CPU backend is auto-detected.
# If the user specifies only a CPU backend, no GPU is used.
# If the user specifies both, we use exactly what was requested.
# ##############################################################################

include(CheckLanguage)
message("")

# Try to compile a minimal CUDA program to verify the compiler actually works.
macro(_try_cuda result)
  try_compile(${result}
    ${CMAKE_BINARY_DIR}/_cuda_check
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/checks/cuda
    cuda_check
    CMAKE_FLAGS "-DCMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}"
  )
endmacro()

# ==============================================================================
# Phase 0: Classify user intent
# ==============================================================================

if(CUDA OR HIP)
  set(_USER_SPECIFIED_GPU TRUE)
else()
  set(_USER_SPECIFIED_GPU FALSE)
endif()

if(OpenMP OR Threads OR Serial)
  set(_USER_SPECIFIED_CPU TRUE)
else()
  set(_USER_SPECIFIED_CPU FALSE)
endif()

if(NOT _USER_SPECIFIED_GPU AND NOT _USER_SPECIFIED_CPU)
  message(
    STATUS
      "${BoldCyan}━━━━━━━━━━━━━━━━━━━━━━━━━━ No device specified, trying to auto-detect ━━━━━━━━━━━━━━━━━━━━━━━━━━━${ColorReset}"
  )
endif()

# ==============================================================================
# Phase 1: Resolve GPU backend (CUDA > HIP > none)
# ==============================================================================

if(_USER_SPECIFIED_GPU)
  # User explicitly requested a GPU backend — validate it.

  if(CUDA)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
      message(
        FATAL_ERROR "CUDA was explicitly requested but no CUDA compiler was found.")
    endif()
    #_try_cuda(_CUDA_COMPILES)
    #if(NOT _CUDA_COMPILES)
    #  message(
    #    FATAL_ERROR "CUDA compiler was found but failed to compile a test program.")
    #endif()
    set(HIP OFF)

  elseif(HIP)
    # Kokkos + HIP: clear offload-arch set by HIP's cmake to avoid conflicts.
    set(GPU_BUILD_TARGETS
        ""
        CACHE STRING "" FORCE)
    set(GPU_TARGETS
        ""
        CACHE STRING "" FORCE)

    check_language(HIP)
    if(NOT CMAKE_HIP_COMPILER)
      message(
        FATAL_ERROR "HIP was explicitly requested but no HIP compiler was found.")
    endif()
    set(CUDA OFF)
    set(HIP_PLATFORM
        "amd"
        CACHE STRING "Set the HIP platform (amd or nvidia)")
  endif()

elseif(NOT _USER_SPECIFIED_CPU)
  # Nothing specified at all — auto-detect GPU.

  # Try CUDA: compiler must exist AND compile a test program.
  check_language(CUDA)
  #if(CMAKE_CUDA_COMPILER)
  #  _try_cuda(_CUDA_COMPILES)
  #endif()

  if(CMAKE_CUDA_COMPILER AND _CUDA_COMPILES)
    set(CUDA ON)
    set(HIP OFF)
  else()
    set(CUDA OFF)

    # Kokkos + HIP: clear offload-arch set by HIP's cmake to avoid conflicts.
    set(GPU_BUILD_TARGETS
        ""
        CACHE STRING "" FORCE)
    set(GPU_TARGETS
        ""
        CACHE STRING "" FORCE)

    check_language(HIP)
    if(CMAKE_HIP_COMPILER)
      set(HIP ON)
      set(HIP_PLATFORM
          "amd"
          CACHE STRING "Set the HIP platform (amd or nvidia)")
    else()
      set(HIP OFF)
    endif()
  endif()

else()
  # User specified only a CPU backend — no GPU.
  set(CUDA OFF)
  set(HIP OFF)
endif()

# ==============================================================================
# Phase 2: Resolve CPU backend (OpenMP > Threads > Serial)
# ==============================================================================

if(_USER_SPECIFIED_CPU)
  # User explicitly requested a CPU backend — validate it.

  if(OpenMP)
    find_package(OpenMP QUIET)
    if(NOT OpenMP_CXX_FOUND)
      message(
        FATAL_ERROR "OpenMP was explicitly requested but was not found.")
    endif()
    set(Threads OFF)
    set(Serial OFF)

  elseif(Threads)
    set(OpenMP OFF)
    set(Serial OFF)

  elseif(Serial)
    set(OpenMP OFF)
    set(Threads OFF)
  endif()

else()
  # Auto-detect CPU backend.

  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    set(OpenMP ON)
    set(Threads OFF)
    set(Serial OFF)
  else()
    set(OpenMP OFF)
    set(Threads ON)
    set(Serial OFF)
  endif()
endif()

# ==============================================================================
# Phase 3: Final validation and summary
# ==============================================================================

if(NOT CUDA
   AND NOT HIP
   AND NOT OpenMP
   AND NOT Threads
   AND NOT Serial)
  message(
    FATAL_ERROR
      "No valid device configuration found. Please specify at least one of: CUDA, HIP, OpenMP, Threads, or Serial."
  )
endif()

set(_GPU_BACKEND "none")
if(CUDA)
  set(_GPU_BACKEND "CUDA")
elseif(HIP)
  set(_GPU_BACKEND "HIP (platform: ${HIP_PLATFORM})")
endif()

set(_CPU_BACKEND "unknown")
if(OpenMP)
  set(_CPU_BACKEND "OpenMP")
elseif(Threads)
  set(_CPU_BACKEND "Threads")
elseif(Serial)
  set(_CPU_BACKEND "Serial")
endif()

message(
  STATUS
    "Device configuration: \n    GPU backend: ${_GPU_BACKEND} \n    CPU backend: ${_CPU_BACKEND}"
)
