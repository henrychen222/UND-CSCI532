cmake_minimum_required (VERSION 2.6)
project (OPENCL_MATRIX_OP)

SET (CMAKE_CXX_FLAGS "-std=c++11 -Wall -O3 -funroll-loops -msse3")

SET (SRC_LIST matrix_op.cpp)

find_package(OpenCL)

message(STATUS "opencl libraries: ${OpenCL_LIBRARIES}")
message(STATUS "opencl include_dirs: ${OpenCL_INCLUDE_DIRS}")

include_directories(${OpenCL_INCLUDE_DIRS})

add_executable(matrix_op matrix_op.cpp)
target_link_libraries(matrix_op ${OpenCL_LIBRARIES})

