message("") # Blank line for better readability
message(
  STATUS
    "${BoldYellow}----------------------- Fetching HDF5 ------------------------${ColorReset}"
)

# enable FetchContent
include(FetchContent)
set(FETCHCONTENT_QUIET
    ON
    CACHE BOOL "")
set(HDF5_ENABLE_PARALLEL
    ${MPI}
    CACHE BOOL "")
set(HDF5_BUILD_EXAMPLES
    OFF
    CACHE BOOL "")
set(HDF5_BUILD_TOOLS
    OFF
    CACHE BOOL "")
set(HDF5_BUILD_HL_LIB
    OFF
    CACHE BOOL "")
set(BUILD_SHARED_LIBS
    OFF
    CACHE BOOL "")
set(BUILD_TESTING
    OFF
    CACHE BOOL "")
# declare HDF5
FetchContent_Declare(
  HDF5 # name of the content
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  URL https://github.com/HDFGroup/hdf5/releases/download/2.0.0/hdf5.tar.gz)
# make available
FetchContent_MakeAvailable(HDF5)

set(HDF5_LIBRARIES
    hdf5-static
    CACHE STRING "")

message(STATUS "HDF5 available")

message(
  STATUS
    "${BoldYellow}--------------------- Fetching HDF5 DONE ---------------------${ColorReset}\n"
)
