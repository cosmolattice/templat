# To find HDF5, we first do a "dry run", i.e. create a new cmake project and try
# to find HDF5 there.
execute_process(
  COMMAND ${CMAKE_COMMAND} -S ${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5/find
          -B ${CMAKE_CURRENT_BINARY_DIR}/cmake/libs/hdf5/find -DMPI=${MPI}
  RESULT_VARIABLE find_result
  OUTPUT_QUIET ERROR_QUIET)

# Now do it for real.
if(find_result EQUAL 0)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libs/hdf5/find/CMakeLists.txt)
endif()
