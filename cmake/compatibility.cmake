if(DEFINED MYFFTW3_PATH)
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${MYFFTW3_PATH}")
endif()

if(DEFINED MYHDF5_PATH)
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${MYHDF5_PATH}")
endif()

if(DEFINED MYFFTW3_PATH OR DEFINED MYHDF5_PATH)
  message(
    WARNING
      "MYFFTW3_PATH and MYHDF5_PATH are deprecated.
  TempLat now automatically downloads and configures FFTW and HDF5 if they are not found on the system, so these variables should no longer be needed.
  If you still wish to use a custom installation of FFTW or HDF5, please set CMAKE_PREFIX_PATH directly."
  )
endif()

if(DEFINED MYPFFT_PATH OR PFFT)
  message(
    FATAL_ERROR
      "PFFT support is no longer available in TempLat. Please use ParaFaFT instead for parallel FFTs (use -DPARAFAFT=ON instead)."
  )
endif()
