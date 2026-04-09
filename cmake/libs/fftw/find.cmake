find_path(FFTW_INCLUDES fftw3.h HINTS ${FFTW_INCLUDES})

# Find the FFTW libraries, threads support is mandatory for TempLat
find_library(FFTW_LIB fftw3 HINTS ${FFTW_LIB_DIR})
find_library(FFTW_THREADS_LIB fftw3_threads HINTS ${FFTW_LIB_DIR})
mark_as_advanced(FFTW_LIB FFTW_THREADS_LIB)

# If float support is enabled, find the single-precision FFTW libraries.
if(FLOAT)
  find_library(FFTWF_LIB fftw3f HINTS ${FFTW_LIB_DIR})
  find_library(FFTWF_THREADS_LIB fftw3f_threads HINTS ${FFTW_LIB_DIR})
  mark_as_advanced(FFTWF_LIB FFTWF_THREADS_LIB)
endif()

# If MPI support is enabled, we need to find the MPI-enabled FFTW libraries.
if(MPI)
  find_library(FFTW_MPI_LIB fftw3_mpi HINTS ${FFTW_LIB_DIR})
  mark_as_advanced(FFTW_MPI_LIB)

  if(FLOAT)
    find_library(FFTWF_MPI_LIB fftw3f_mpi HINTS ${FFTW_LIB_DIR})
    mark_as_advanced(FFTWF_MPI_LIB)
  endif()
endif()

# Catch an edge case: An enabled OpenMP backend without an OpenMP-enabled FFTW.
# In this case, we must also find the OpenMP-enabled FFTW libraries.
if(OPENMP)
  find_library(FFTW_OPENMP_LIB fftw3_omp HINTS ${FFTW_LIB_DIR})
  mark_as_advanced(FFTW_OPENMP_LIB)

  if(FLOAT)
    find_library(FFTWF_OPENMP_LIB fftw3f_omp HINTS ${FFTW_LIB_DIR})
    mark_as_advanced(FFTWF_OPENMP_LIB)
  endif()
endif()

# ##############################################################################
# Assemble FFTW_LIBRARIES in correct dependency order (dependents first). With
# static archives the linker resolves symbols left-to-right, so libraries that
# *use* symbols must appear before the libraries that *provide* them: fftw3_mpi
# / fftw3_omp  →  fftw3_threads  →  fftw3  (and likewise for float)
# ##############################################################################

set(FFTW_LIBRARIES "")

# MPI wrappers depend on the core libraries, so they must come first.
if(MPI)
  list(APPEND FFTW_LIBRARIES ${FFTW_MPI_LIB})
  if(FLOAT)
    list(APPEND FFTW_LIBRARIES ${FFTWF_MPI_LIB})
  endif()
endif()

# OpenMP wrappers also depend on the core libraries.
if(OPENMP)
  list(APPEND FFTW_LIBRARIES ${FFTW_OPENMP_LIB})
  if(FLOAT)
    list(APPEND FFTW_LIBRARIES ${FFTWF_OPENMP_LIB})
  endif()
endif()

# Threads and core come last (threads depends on core).
list(APPEND FFTW_LIBRARIES ${FFTW_THREADS_LIB} ${FFTW_LIB})
if(FLOAT)
  list(APPEND FFTW_LIBRARIES ${FFTWF_THREADS_LIB} ${FFTWF_LIB})
endif()

# And finally, check that we found everything we needed
set(_FFTW_REQUIRED_VARS FFTW_INCLUDES FFTW_LIB FFTW_THREADS_LIB)
if(FLOAT)
  list(APPEND _FFTW_REQUIRED_VARS FFTWF_LIB FFTWF_THREADS_LIB)
endif()
if(MPI)
  list(APPEND _FFTW_REQUIRED_VARS FFTW_MPI_LIB)
  if(FLOAT)
    list(APPEND _FFTW_REQUIRED_VARS FFTWF_MPI_LIB)
  endif()
endif()
if(OPENMP)
  list(APPEND _FFTW_REQUIRED_VARS FFTW_OPENMP_LIB)
  if(FLOAT)
    list(APPEND _FFTW_REQUIRED_VARS FFTWF_OPENMP_LIB)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG ${_FFTW_REQUIRED_VARS})

# extract FFTW_LIB_DIR from FFTW_LIB
if(FFTW_LIB)
  get_filename_component(FFTW_LIB_DIR "${FFTW_LIB}" DIRECTORY)
endif()
