# Glob all test sources and then create targets for each of them
file(GLOB bench_sources "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
foreach(fn ${bench_sources})
  get_filename_component(f_name ${fn} NAME) # extract the filename from the full
                                            # path

  # remove cpp and create a string appending test_ at the begining. Allows to
  # create a rule test_filename
  string(REPLACE ".cpp" "" target_name "bench-${f_name}")

  add_executable(${target_name} ${fn})
  target_link_libraries(${target_name} TempLat::TempLat)

  # Grid/step overrides for the halo_scaling benchmark (compile-time constants).
  # Set e.g. -DHALO_NGRID=256 -DHALO_NSTEPS=100 at configure time.
  if(DEFINED HALO_NGRID)
    target_compile_definitions(${target_name} PRIVATE HALO_NGRID=${HALO_NGRID})
  endif()
  if(DEFINED HALO_NSTEPS)
    target_compile_definitions(${target_name} PRIVATE HALO_NSTEPS=${HALO_NSTEPS})
  endif()
endforeach()

# ##############################################################################
# Codegen guardrail: the CPU lattice kernels must actually vectorize.
#
# Nothing else in the suite can see this. A scalar kernel is correct, just ~2.5x slow, so it passes
# every test we have -- which is how a 100%-scalar SU(2) kick survived. See PERFORMANCE.md.
#
# The vectorized dispatch is opt-in (PERFORMANCE.md, section 5), so this test needs its own build with
# it enabled: against the default binary it would find no inner-loop kernels and fail for the wrong
# reason. It still guards the path, so that it can be evaluated again when it is picked up.
# ##############################################################################

if(TEMPLAT_TEST
   AND NOT CUDA
   AND NOT HIP
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
  add_executable(bench-su2_evolution_vec "${CMAKE_CURRENT_SOURCE_DIR}/su2_evolution.cpp")
  target_link_libraries(bench-su2_evolution_vec TempLat::TempLat)
  target_compile_definitions(bench-su2_evolution_vec PRIVATE TEMPLAT_VECTORIZED_MDRANGE)

  add_test(
    NAME test-vectorization
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_vectorization.sh --cpu-build
            $<TARGET_FILE:bench-su2_evolution_vec> 100)
  set_tests_properties(test-vectorization PROPERTIES LABELS "benchmarks;codegen")
endif()
