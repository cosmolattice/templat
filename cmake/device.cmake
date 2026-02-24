# ##############################################################################
# Find out which backends are available
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
      "${Yellow}---------- No device specified, trying to auto-detect ----------${ColorReset}"
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
   AND SERIAL)
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

# ##############################################################################
# Choose a device provider
# ##############################################################################

if(DEVICE_PROVIDER STREQUAL "Kokkos")
  set(KOKKOS ON)

  message(
    STATUS
      "${Green}---------- Using Kokkos as device provider ----------${ColorReset}\n"
  )

  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/device/kokkos.cmake)
  if(KOKKOSFFT)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/device/kokkos-fft.cmake)
  endif()

  set(CMAKE_REQUIRED_QUIET ON)
  set(CMAKE_MESSAGE_LOG_LEVEL WARNING)
  set(CMAKE_MESSAGE_LOG_LEVEL STATUS)
else()
  message(FATAL_ERROR "Unknown DEVICE_PROVIDER option: ${DEVICE_PROVIDER}.
      Supported options: Kokkos")
endif()

function(target_link_device target)
  if(NOT TARGET ${target})
    message(FATAL_ERROR "Target ${target} does not exist.")
  endif()

  # Link against the device provider's libraries and define a preprocessor macro
  # to indicate which device provider is being used.
  if(KOKKOS)
    target_link_libraries(${target} INTERFACE Kokkos::kokkos)
    target_compile_definitions(${target} INTERFACE DEVICE_KOKKOS)
  endif()

  # We handle here also KokkosFFT
  if(KOKKOSFFT)
    target_link_libraries(${target} INTERFACE KokkosFFT::fft)
    target_compile_definitions(${target} INTERFACE HAVE_KOKKOSFFT)
  endif()

  # Define a preprocessor macro to indicate which device is being used. This can
  # be used in the code to conditionally compile device-specific code.
  if(CUDA)
    target_compile_definitions(TempLat INTERFACE DEVICE_CUDA)
  elseif(HIP)
    target_compile_definitions(TempLat INTERFACE DEVICE_HIP)
  else()
    target_compile_definitions(TempLat INTERFACE DEVICE_CPU)
  endif()
endfunction()
