# ##############################################################################
# FLAGS
# ##############################################################################

include(ProcessorCount)
ProcessorCount(CMAKE_PROCESSOR_COUNT)

set(FFTW_CONFIGFLAGS
    " --enable-threads --disable-shared --enable-static --enable-sse2 --enable-avx --enable-avx2 --enable-avx512"
    CACHE
      STRING
      "Configuration flags for FFTW. Default is --enable-threads --disable-shared --enable-static"
)

if(MPI)
  set(FFTW_CONFIGFLAGS "${FFTW_CONFIGFLAGS} --enable-mpi")
endif()

if(OPENMP)
  set(FFTW_CONFIGFLAGS "${FFTW_CONFIGFLAGS} --enable-openmp")
endif()

# ##############################################################################
# Fetch FFTW tarball and extract
# ##############################################################################

message("") # Blank line for better readability
message(
  "${BoldYellow}┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching FFTW ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓${ColorReset}"
)

# Download FFTW
set(FFTW_URL "https://www.fftw.org/fftw-3.3.10.tar.gz")
set(FFTW_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/fftw")

if(NOT EXISTS "${FFTW_DIR}/fftw-3.3.10.tar.gz")
  message(STATUS "Downloading FFTW from ${FFTW_URL}...")
  file(DOWNLOAD "${FFTW_URL}" "${FFTW_DIR}/fftw-3.3.10.tar.gz")
else()
  message(
    STATUS
      "FFTW archive already exists at ${FFTW_DIR}/fftw-3.3.10.tar.gz, skipping download."
  )
endif()

# Extract FFTW
if(NOT EXISTS "${FFTW_DIR}/fftw-3.3.10")
  message(STATUS "Extracting FFTW...")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xzf "${FFTW_DIR}/fftw-3.3.10.tar.gz"
    WORKING_DIRECTORY "${FFTW_DIR}"
    RESULT_VARIABLE fftw_configure_ret)
  if(NOT fftw_configure_ret EQUAL 0)
    message(FATAL_ERROR "FFTW extraction failed.")
  endif()
else()
  message(
    STATUS
      "FFTW already extracted at ${FFTW_DIR}/fftw-3.10, skipping extraction.")
endif()

# ##############################################################################
# Build FFTW (double)
# ##############################################################################

message(STATUS "Configuring FFTW (double)... (flags: ${FFTW_CONFIGFLAGS})")
execute_process(
  COMMAND bash -c "./configure --prefix=${FFTW_DIR} ${FFTW_CONFIGFLAGS}"
  WORKING_DIRECTORY "${FFTW_DIR}/fftw-3.3.10/"
  OUTPUT_FILE "${FFTW_DIR}/configure.log"
  ERROR_FILE "${FFTW_DIR}/configure.log"
  RESULT_VARIABLE fftw_configure_ret)
if(NOT fftw_configure_ret EQUAL 0)
  message(
    FATAL_ERROR
      "FFTW configure failed. See ${FFTW_DIR}/configure.log for details.")
endif()

message(
  STATUS "Building FFTW (double) with ${CMAKE_PROCESSOR_COUNT} threads...")
execute_process(
  COMMAND make -j ${CMAKE_PROCESSOR_COUNT}
  WORKING_DIRECTORY "${FFTW_DIR}/fftw-3.3.10/"
  OUTPUT_FILE "${FFTW_DIR}/build.log"
  ERROR_FILE "${FFTW_DIR}/build.log"
  RESULT_VARIABLE fftw_build_ret)
if(NOT fftw_build_ret EQUAL 0)
  message(
    FATAL_ERROR "FFTW build failed. See ${FFTW_DIR}/build.log for details.")
endif()

message(STATUS "Installing FFTW (double)...")
execute_process(
  COMMAND make install
  WORKING_DIRECTORY "${FFTW_DIR}/fftw-3.3.10/"
  OUTPUT_FILE "${FFTW_DIR}/install.log"
  ERROR_FILE "${FFTW_DIR}/install.log"
  RESULT_VARIABLE fftw_install_ret)
if(NOT fftw_install_ret EQUAL 0)
  message(
    FATAL_ERROR "FFTW install failed. See ${FFTW_DIR}/install.log for details.")
endif()

# ##############################################################################
# Build FFTW (float)
# ##############################################################################

if(FLOAT)
  message(
    STATUS
      "Configuring FFTW (float)... (flags: --enable-float --enable-sse ${FFTW_CONFIGFLAGS})"
  )
  execute_process(
    COMMAND
      bash -c
      "./configure --prefix=${FFTW_DIR} --enable-float --enable-sse ${FFTW_CONFIGFLAGS}"
    WORKING_DIRECTORY "${FFTW_DIR}/fftw-3.3.10/"
    OUTPUT_FILE "${FFTW_DIR}/configure.log"
    ERROR_FILE "${FFTW_DIR}/configure.log"
    RESULT_VARIABLE fftw_configure_ret)
  if(NOT fftw_configure_ret EQUAL 0)
    message(
      FATAL_ERROR
        "FFTWf configure failed. See ${FFTW_DIR}/configure.log for details.")
  endif()

  message(
    STATUS "Building FFTW (float) with ${CMAKE_PROCESSOR_COUNT} threads...")
  execute_process(
    COMMAND make -j ${CMAKE_PROCESSOR_COUNT}
    WORKING_DIRECTORY "${FFTW_DIR}/fftw-3.3.10/"
    OUTPUT_FILE "${FFTW_DIR}/build.log"
    ERROR_FILE "${FFTW_DIR}/build.log"
    RESULT_VARIABLE fftw_build_ret)
  if(NOT fftw_build_ret EQUAL 0)
    message(
      FATAL_ERROR "FFTWf build failed. See ${FFTW_DIR}/build.log for details.")
  endif()

  message(STATUS "Installing FFTW (float)...")
  execute_process(
    COMMAND make install
    WORKING_DIRECTORY "${FFTW_DIR}/fftw-3.3.10/"
    OUTPUT_FILE "${FFTW_DIR}/install.log"
    ERROR_FILE "${FFTW_DIR}/install.log"
    RESULT_VARIABLE fftw_install_ret)
  if(NOT fftw_install_ret EQUAL 0)
    message(
      FATAL_ERROR
        "FFTWf install failed. See ${FFTW_DIR}/install.log for details.")
  endif()
endif()

# ##############################################################################
# Finish
# ##############################################################################

set(CMAKE_PREFIX_PATH "${FFTW_DIR};${CMAKE_PREFIX_PATH}")

set(FFTW_INCLUDES
    "${FFTW_DIR}/include"
    CACHE PATH "Path to FFTW include directory" FORCE)
set(FFTW_LIB_DIR
    "${FFTW_DIR}/lib"
    CACHE PATH "Path to FFTW library directory" FORCE)

message(
  "${BoldYellow}┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching FFTW DONE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛${ColorReset}\n"
)
