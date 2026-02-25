include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5/find.cmake)

if(NOT HDF5_FOUND OR NOT CORRECT_HDF_FOUND)
  if(${AUTOBUILD_HDF5})
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5/get.cmake)
  else()
    message(
      FATAL_ERROR
        "HDF5 libraries not found. Set AUTOBUILD_HDF5 to ON to automatically build from source."
    )
  endif()
endif()
