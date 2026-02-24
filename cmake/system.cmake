# ##############################################################################
# Setting compiler flags depending on the system, build type and compiler.
# ##############################################################################

set(NATIVE
    (NOT APPLE)
    CACHE STRING "Set whether to use -march=native.")
set(SSE
    ON
    CACHE
      BOOL
      "Set the sse instruction flag for the compiler. Default is OFF, but most likely you want to set is to ON."
)
set(AVX
    OFF
    CACHE
      STRING
      "Set the avx instruction flags during compilation. Can choose between OFF, mavx, mavx2, mavx512f (default = OFF). If you want something else, feel free to change the CMakeLists."
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(NATIVE)
    target_compile_options(TempLat INTERFACE -march=native)
  endif()
  target_compile_options(TempLat INTERFACE -ffast-math)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(TempLat INTERFACE -Wno-nan-infinity-disabled)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(TempLat INTERFACE -fsanitize=address)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(TempLat INTERFACE -fconcepts-diagnostics-depth=8)
  endif()
endif()

if(${SSE})
  target_compile_options(TempLat INTERFACE -msse -msse2 -msse3 -msse4)
endif()

if(${AVX} STREQUAL mavx)
  target_compile_options(TempLat INTERFACE -mavx)
elseif(${AVX} STREQUAL mavx2)
  target_compile_options(TempLat INTERFACE -mavx -mavx2)
elseif(${AVX} STREQUAL mavx512f)
  target_compile_options(TempLat INTERFACE -mavx -mavx2 -mavx512f)
elseif(NOT ${AVX} STREQUAL OFF)
  message(FATAL_ERROR "AVX=${AVX} is not a valid option.")
endif()

# Make sure GCC can handle huge amounts of constexpr computations.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(TempLat INTERFACE -fconstexpr-loop-limit=2621440
                                           -fconstexpr-ops-limit=335544320)
endif()

target_compile_options(TempLat INTERFACE -Wall -Wshadow)

# ##############################################################################
# These flags are deprecated from last version, notify the user if they are
# still using them.
# ##############################################################################

if(DEFINED G++OPT)
  message(
    FATAL_ERROR
      "G++OPT is deprecated. Please use CMAKE_BUILD_TYPE=Release or Debug instead, and set NATIVE, SSE, and AVX as needed."
  )
endif()
if(DEFINED G++AVX)
  message(
    FATAL_ERROR
      "G++AVX is deprecated. Please use CMAKE_BUILD_TYPE=Release or Debug instead, and set NATIVE, SSE, and AVX as needed."
  )
endif()
if(DEFINED G++SSE)
  message(
    FATAL_ERROR
      "G++SSE is deprecated. Please use CMAKE_BUILD_TYPE=Release or Debug instead, and set NATIVE, SSE, and AVX as needed."
  )
endif()
