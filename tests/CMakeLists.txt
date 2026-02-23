# ---------- Number of MPI processes ----------
set(NPROCESSES
    4
    CACHE STRING "Number of MPI processes used for testing")

# ---------- Helper function ----------
function(add_templat_test execName mainFile)
  add_executable(${execName} ${mainFile})
  target_link_libraries(${execName} PUBLIC TempLat::TempLat)
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
    add_test(
      NAME ${target_name}
      COMMAND
        bash -c
        "make ${target_name} &&
                ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROCESSES} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${target_name}> ${MPIEXEC_POSTFLAGS}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  else()
    add_test(
      NAME ${target_name}
      COMMAND bash -c "make ${target_name} && $<TARGET_FILE:${target_name}>"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  set_tests_properties(${target_name} PROPERTIES LABELS "${folder_labels}")
endforeach()
