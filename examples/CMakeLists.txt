cmake_minimum_required(VERSION 3.16)
project(my_project LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
  TempLat
  GIT_REPOSITORY https://github.com/cosmolattice/templat.git
  GIT_TAG main)
FetchContent_MakeAvailable(TempLat)

# Add executable and link it against TempLat
add_executable(langevin langevin.cpp)
target_link_libraries(langevin PRIVATE TempLat::TempLat)

add_executable(average average.cpp)
target_link_libraries(average PRIVATE TempLat::TempLat)

add_executable(classical_phi4 classical_phi4.cpp)
target_link_libraries(classical_phi4 PRIVATE TempLat::TempLat)

add_executable(classical_su2 classical_su2.cpp)
target_link_libraries(classical_su2 PRIVATE TempLat::TempLat)
