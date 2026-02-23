# TempLat Test Build System
# This mirrors cosmolattice's src/tests/CMakeLists.txt but is self-contained.

# ---------- Dependencies ----------
include(${CMAKE_SOURCE_DIR}/cmake/device.cmake)

find_package(Threads REQUIRED)

option(MPI "Enable MPI support" OFF)
if(MPI)
    find_package(MPI COMPONENTS CXX REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_MPI")
endif()

set(MYFFTW3_PATH "" CACHE STRING "Path to FFTW3 installation")
include(${CMAKE_SOURCE_DIR}/cmake/libs/fftw.cmake)

option(HDF5 "Enable HDF5 support" OFF)
if(HDF5)
    find_package(HDF5 REQUIRED)
    include_directories(${HDF5_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_HDF5")
endif()

option(PARAFAFT "Enable ParaFaFT support" OFF)
if(PARAFAFT AND MPI)
    include(${CMAKE_SOURCE_DIR}/cmake/libs/parafaft.cmake)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_PARAFAFT")
endif()

if(NOT DEFINED FFTW_LIB AND NOT DEFINED KOKKOSFFT AND NOT DEFINED PARAFAFT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOFFT")
endif()

# ---------- Number of MPI processes ----------
set(NPROCESSES 4 CACHE STRING "Number of MPI processes used for testing")

# ---------- Helper function ----------
function(add_templat_test execName mainFile)
    add_executable(${execName} ${mainFile})
    target_link_libraries(${execName} PUBLIC TempLat::TempLat)
    target_link_libraries(${execName} PUBLIC Threads::Threads)
    target_link_libraries(${execName} PUBLIC ${FFTW_LIBRARIES})
    target_link_device(${execName})
    if(MPI)
        target_link_libraries(${execName} PUBLIC ${MPI_CXX_LIBRARIES})
    endif()
    if(PARAFAFT AND MPI)
        target_link_libraries(${execName} PUBLIC parafaft::parafaft)
    endif()
    if(HDF5)
        target_link_libraries(${execName} PUBLIC ${HDF5_LIBRARIES})
    endif()
endfunction()

# ---------- Discover and register tests ----------
file(GLOB_RECURSE test_sources "${CMAKE_CURRENT_SOURCE_DIR}/TempLat/*/*.cpp"
                                "${CMAKE_CURRENT_SOURCE_DIR}/TempLat/*.cpp")

foreach(fn ${test_sources})
    get_filename_component(f_name ${fn} NAME)
    string(REPLACE ".cpp" "" test_name "${f_name}")
    set(target_name "test-${test_name}")

    add_templat_test(${target_name} ${CMAKE_CURRENT_SOURCE_DIR}/tddmain.cpp)
    target_sources(${target_name} PRIVATE ${fn})
    target_compile_definitions(${target_name} PUBLIC CHECKBOUNDS)

    # Labels for ctest filtering
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" relative_path "${fn}")
    get_filename_component(folder_path "${relative_path}" DIRECTORY)
    string(REPLACE "/" ";" folder_labels "${folder_path}")

    if(MPI)
        add_test(NAME ${target_name}
            COMMAND bash -c "make ${target_name} &&
                ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROCESSES} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${target_name}> ${MPIEXEC_POSTFLAGS}"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    else()
        add_test(NAME ${target_name}
            COMMAND bash -c "make ${target_name} && $<TARGET_FILE:${target_name}>"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    set_tests_properties(${target_name} PROPERTIES LABELS "${folder_labels}")
endforeach()
