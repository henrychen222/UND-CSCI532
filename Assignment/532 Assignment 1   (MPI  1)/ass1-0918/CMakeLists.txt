cmake_minimum_required (VERSION 2.6)
project (ass1)

SET (CMAKE_C_FLAGS "-Wall -std=c99")

add_executable(kmeans_non_mpi kmeans_non_mpi.c)
target_link_libraries(kmeans_non_mpi m)
