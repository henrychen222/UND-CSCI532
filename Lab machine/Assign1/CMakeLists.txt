cmake_minimum_required (VERSION 2.6)
project (ass1)

# The version number.
(ass1_VERSION_MAJOR 1)
(ass1_VERSION_MINOR 0)

(CMAKE_C_COMPILER	"/usr/bin/clang")
(CMAKE_C_FLAGS	"-Wall -std=c99")

(CMAKE_CXX_COMPILER	"/usr/bin/clang++")
(CMAKE_CXX_FLAGS	"-Wall")

(CMAKE_AR      "/usr/bin/llvm-ar")
(CMAKE_LINKER  "/usr/bin/llvm-ld")
(CMAKE_NM      "/usr/bin/llvm-nm")
(CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
(CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

add_executable(kmeans-mpi kmeans.c)
target_link_libraries(kmeans-mpi m)
