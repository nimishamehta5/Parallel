cmake_minimum_required(VERSION 3.9.1)
project(p3)
find_package(CUDA QUIET REQUIRED)
find_package(MPI REQUIRED)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Werror")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
#set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Werror")

include_directories( ${MPI_CXX_INCLUDE_PATH})
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})

include_directories (cuda_codes,mpi_codes,thread_codes)

set( SOURCE_FILES thread_codes/thread_final.cc)
#thread_codes/complex.h thread_codes/complex.cc thread_codes/input_image.h thread_codes/input_image.cc
set(SOURCE_FILES2 mpi_codes/mpi_final.cc)

set(SOURCE_FILES3 cuda_codes/cuda_final.cu)

add_executable(p31 ${SOURCE_FILES})
cuda_add_executable(p33 ${SOURCE_FILES3})
add_executable( p32 ${SOURCE_FILES2})
target_link_libraries(p32 ${MPI_CXX_LIBRARIES})
