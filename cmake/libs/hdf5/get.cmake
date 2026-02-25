message("") # Blank line for better readability
message(
  STATUS
    "${BoldYellow}┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching HDF5 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓${ColorReset}"
)

# enable FetchContent
include(FetchContent)
set(FETCHCONTENT_QUIET
    ON
    CACHE BOOL "")
mark_as_advanced(FETCHCONTENT_QUIET)
set(HDF5_ENABLE_PARALLEL
    ${MPI}
    CACHE BOOL "")
set(HDF5_BUILD_EXAMPLES
    OFF
    CACHE BOOL "")
mark_as_advanced(HDF5_BUILD_EXAMPLES)
set(HDF5_BUILD_TOOLS
    OFF
    CACHE BOOL "")
mark_as_advanced(HDF5_BUILD_TOOLS)
set(HDF5_BUILD_HL_LIB
    OFF
    CACHE BOOL "")
mark_as_advanced(HDF5_BUILD_HL_LIB)
set(BUILD_SHARED_LIBS
    OFF
    CACHE BOOL "")
mark_as_advanced(BUILD_SHARED_LIBS)
set(BUILD_TESTING
    OFF
    CACHE BOOL "")
mark_as_advanced(BUILD_TESTING)
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
    "${BoldYellow}┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching HDF5 DONE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛${ColorReset}\n"
)
