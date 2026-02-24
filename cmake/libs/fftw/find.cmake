find_path(FFTW_INCLUDES fftw3.h)
mark_as_advanced(FFTW_INCLUDES)

# Find the FFTW libraries, threads support is mandatory for TempLat
find_library(FFTW_LIB fftw3)
find_library(FFTW_THREADS_LIB fftw3_threads)
mark_as_advanced(FFTW_LIB FFTW_THREADS_LIB)

# If float support is enabled, find the single-precision FFTW libraries.
if(Float)
  find_library(FFTWF_LIB fftw3f)
  find_library(FFTWF_THREADS_LIB fftw3f_threads)
  mark_as_advanced(FFTWF_LIB FFTWF_THREADS_LIB)
endif()

# If MPI support is enabled, we need to find the MPI-enabled FFTW libraries.
if(MPI)
  find_library(FFTW_MPI_LIB fftw3_mpi)
  mark_as_advanced(FFTW_MPI_LIB)

  if(Float)
    find_library(FFTWF_MPI_LIB fftw3f_mpi)
    mark_as_advanced(FFTWF_MPI_LIB)
  endif()
endif()

# Catch an edge case: An enabled OpenMP backend without an OpenMP-enabled FFTW.
# In this case, we must also find the OpenMP-enabled FFTW libraries.
if(OpenMP)
  find_library(FFTW_OPENMP_LIB fftw3_omp)
  mark_as_advanced(FFTW_OPENMP_LIB)

  if(Float)
    find_library(FFTWF_OPENMP_LIB fftw3f_omp)
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
  if(Float)
    list(APPEND FFTW_LIBRARIES ${FFTWF_MPI_LIB})
  endif()
endif()

# OpenMP wrappers also depend on the core libraries.
if(OpenMP)
  list(APPEND FFTW_LIBRARIES ${FFTW_OPENMP_LIB})
  if(Float)
    list(APPEND FFTW_LIBRARIES ${FFTWF_OPENMP_LIB})
  endif()
endif()

# Threads and core come last (threads depends on core).
list(APPEND FFTW_LIBRARIES ${FFTW_THREADS_LIB} ${FFTW_LIB})
if(Float)
  list(APPEND FFTW_LIBRARIES ${FFTWF_THREADS_LIB} ${FFTWF_LIB})
endif()

# And finally, check that we found everything we needed
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDES FFTW_LIBRARIES)
