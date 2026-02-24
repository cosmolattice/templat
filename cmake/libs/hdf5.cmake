include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5/find.cmake)

if(NOT HDF5_FOUND OR NOT CORRECT_HDF_FOUND)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5/get.cmake)
endif()
