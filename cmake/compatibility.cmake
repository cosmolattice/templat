# ##############################################################################
# Device option renames (TempLat 1.0 → 1.1)
# ##############################################################################

if(DEFINED OpenMP AND NOT DEFINED OPENMP)
  message(
    WARNING
      "${BoldMagenta}OpenMP is deprecated. Please use OPENMP instead.${ColorReset}"
  )
  set(OPENMP ${OpenMP})
endif()

if(DEFINED Threads AND NOT DEFINED PTHREADS)
  message(
    WARNING
      "${BoldMagenta}Threads is deprecated. Please use PTHREADS instead.${ColorReset}"
  )
  set(PTHREADS ${Threads})
endif()

if(DEFINED Serial AND NOT DEFINED NOTHREADING)
  message(
    WARNING
      "${BoldMagenta}Serial is deprecated. Please use NOTHREADING instead.${ColorReset}"
  )
  set(NOTHREADING ${Serial})
endif()

if(DEFINED Float AND NOT DEFINED FLOAT)
  message(
    WARNING
      "${BoldMagenta}Float is deprecated. Please use FLOAT instead.${ColorReset}"
  )
  set(FLOAT ${Float})
endif()

# ##############################################################################
# Dependency handling from CL1.0
# ##############################################################################

if(DEFINED MYFFTW3_PATH)
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${MYFFTW3_PATH}")
endif()

if(DEFINED MYHDF5_PATH)
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${MYHDF5_PATH}")
endif()

if(DEFINED MYFFTW3_PATH OR DEFINED MYHDF5_PATH)
  message(
    WARNING
      "${BoldMagenta}MYFFTW3_PATH and MYHDF5_PATH are deprecated. TempLat now automatically downloads and configures FFTW and HDF5 if they are not found on the system, so these variables should no longer be needed. If you still wish to use a custom installation of FFTW or HDF5, please set CMAKE_PREFIX_PATH directly.${ColorReset}"
  )
endif()

# PFFT support was removed in CL1.0, so if the user is still trying to use it,
# notify them to switch to ParaFaFT instead.

if(DEFINED MYPFFT_PATH OR PFFT)
  message(
    FATAL_ERROR
      "PFFT support is no longer available in TempLat. Please use ParaFaFT instead for parallel FFTs (use -DPARAFAFT=ON instead)."
  )
endif()

# ##############################################################################
# These flags are deprecated from last version, notify the user if they are
# still using them.
# ##############################################################################

if(DEFINED G++OPT)
  message(
    WARNING
      "${BoldMagenta}G++OPT is deprecated. Please use CMAKE_BUILD_TYPE=Release or Debug instead, and set NATIVE, SSE, and AVX as needed.${ColorReset}"
  )
  if(${G++OPT} STREQUAL "G")
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING "Build type" FORCE)
    message(
      WARNING
        "${BoldMagenta}G++OPT=G translated to setting CMAKE_BUILD_TYPE to Debug, which enables debug flags and disables optimizations.${ColorReset}"
    )
  elseif(${G++OPT} STREQUAL "O1")
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING "Build type" FORCE)
    message(
      WARNING
        "${BoldMagenta}G++OPT=O1 translated to setting CMAKE_BUILD_TYPE to Release, which enables optimizations at O3 level.${ColorReset}"
    )
  elseif(${G++OPT} STREQUAL "O2")
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING "Build type" FORCE)
    message(
      WARNING
        "${BoldMagenta}G++OPT=O2 translated to setting CMAKE_BUILD_TYPE to Release, which enables optimizations at O3 level.${ColorReset}"
    )
  elseif(${G++OPT} STREQUAL "O3")
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING "Build type" FORCE)
    message(
      WARNING
        "${BoldMagenta}G++OPT=O3 translated to setting CMAKE_BUILD_TYPE to Release, which enables optimizations at O3 level.${ColorReset}"
    )
  elseif(${G++OPT} STREQUAL "Ofast")
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING "Build type" FORCE)
    message(
      WARNING
        "${BoldMagenta}G++OPT=Ofast translated to setting CMAKE_BUILD_TYPE to Release, which enables optimizations at O3 level.${ColorReset}"
    )
  else()
    message(
      WARNING
        "${BoldMagenta}G++OPT=${G++OPT} has never been a valid option.${ColorReset}"
    )
  endif()
endif()

if(DEFINED G++AVX)
  message(
    WARNING
      "${BoldMagenta}G++AVX is deprecated. Please use CMAKE_BUILD_TYPE=Release or Debug instead, and set NATIVE, SSE, and AVX as needed.${ColorReset}"
  )
  set(AVX
      ${G++AVX}
      CACHE
        STRING
        "Set the avx instruction flags during compilation. Can choose between OFF, mavx, mavx2, mavx512f (default = OFF)."
        FORCE)
endif()

if(DEFINED G++SSE)
  message(
    WARNING
      "${BoldMagenta}G++SSE is deprecated. Please use CMAKE_BUILD_TYPE=Release or Debug instead, and set NATIVE, SSE, and AVX as needed.${ColorReset}"
  )
  set(SSE
      ON
      CACHE
        BOOL
        "Set the sse instruction flag for the compiler. Default is OFF, but most likely you want to set is to ON."
        FORCE)
endif()
