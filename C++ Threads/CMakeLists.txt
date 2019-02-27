cmake_minimum_required(VERSION 3.8)
project(se3)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11 -pthread")
find_package (Threads)
include_directories(
src
)
set(SOURCE_FILES
src/heat2D3D.cc
)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
