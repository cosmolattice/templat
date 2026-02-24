find_path(FFTW_INCLUDES fftw3.h)
mark_as_advanced(FFTW_INCLUDES)

# Find the FFTW libraries, threads support is mandatory for TempLat
find_library(FFTW_LIB fftw3)
find_library(FFTW_THREADS_LIB fftw3_threads)
set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTW_THREADS_LIB})
mark_as_advanced(FFTW_LIB FFTW_THREADS_LIB)

# If float support is enabled, find the single-precision FFTW libraries.
if(Float)
  find_library(FFTWF_LIB fftw3f)
  find_library(FFTWF_THREADS_LIB fftw3f_threads)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWF_LIB} ${FFTWF_THREADS_LIB})
  mark_as_advanced(FFTWF_LIB FFTWF_THREADS_LIB)
endif()

# If MPI support is enabled, we need to find the MPI-enabled FFTW libraries.
if(MPI)
  find_library(FFTW_MPI_LIB fftw3_mpi)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_MPI_LIB})
  mark_as_advanced(FFTW_MPI_LIB)

  if(Float)
    find_library(FFTWF_MPI_LIB fftw3f_mpi)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWF_MPI_LIB})
    mark_as_advanced(FFTWF_MPI_LIB)
  endif()
endif()

# Catch an edge case: An enabled OpenMP backend without an OpenMP-enabled FFTW.
# In this case, we must also find the OpenMP-enabled FFTW libraries.
if(OpenMP AND KOKKOSFFT)
  find_library(FFTW_OPENMP_LIB fftw3_omp)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_OPENMP_LIB})
  mark_as_advanced(FFTW_OPENMP_LIB)

  if(Float)
    find_library(FFTWF_OPENMP_LIB fftw3f_omp)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWF_OPENMP_LIB})
    mark_as_advanced(FFTWF_OPENMP_LIB)
  endif()
endif()

# And finally, check that we found everything we needed
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDES FFTW_LIBRARIES)
