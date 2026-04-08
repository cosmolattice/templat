cmake_minimum_required(VERSION 3.18)
project(check_gpu_mpi LANGUAGES CXX CUDA)

find_package(MPI REQUIRED)

add_executable(check_gpu_mpi check_gpu_mpi.cu)
target_link_libraries(check_gpu_mpi PRIVATE MPI::MPI_CXX)
target_compile_features(check_gpu_mpi PRIVATE cxx_std_17)
