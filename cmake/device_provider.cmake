# ##############################################################################
# Choose a device provider
# ##############################################################################

if(DEVICE_PROVIDER STREQUAL "Kokkos")
  set(KOKKOS ON)

  message(
    STATUS
      "${BoldGreen}━━━━━━━━━━━━━━━━━━━━━━━━ Using Kokkos as device provider ━━━━━━━━━━━━━━━━━━━━━━━━━${ColorReset}\n"
  )

  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/kokkos.cmake)

  set(CMAKE_REQUIRED_QUIET ON)
  set(CMAKE_MESSAGE_LOG_LEVEL WARNING)
  set(CMAKE_MESSAGE_LOG_LEVEL STATUS)
else()
  message(FATAL_ERROR "Unknown DEVICE_PROVIDER option: ${DEVICE_PROVIDER}.
      Supported options: Kokkos")
endif()

# ##############################################################################
# Link and define preprocessor macros
# ##############################################################################

# Link against the device provider's libraries and define a preprocessor macro
# to indicate which device provider is being used.
if(KOKKOS)
  target_link_libraries(TempLat INTERFACE Kokkos::kokkos)
  target_compile_definitions(TempLat INTERFACE DEVICE_KOKKOS)
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
