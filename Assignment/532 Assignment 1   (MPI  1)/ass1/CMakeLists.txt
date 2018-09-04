cmake_minimum_required (VERSION 2.6)
project (ass1)

# The version number.
set (ass1_VERSION_MAJOR 1)
set (ass1_VERSION_MINOR 0)

SET (CMAKE_C_COMPILER	"/usr/bin/clang")
SET (CMAKE_C_FLAGS	"-Wall -std=c99")

SET (CMAKE_CXX_COMPILER	"/usr/bin/clang++")
SET (CMAKE_CXX_FLAGS	"-Wall")

SET (CMAKE_AR      "/usr/bin/llvm-ar")
SET (CMAKE_LINKER  "/usr/bin/llvm-ld")
SET (CMAKE_NM      "/usr/bin/llvm-nm")
SET (CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
SET (CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

add_executable(kmeans-mpi kmeans.c)
target_link_libraries(kmeans-mpi m)
