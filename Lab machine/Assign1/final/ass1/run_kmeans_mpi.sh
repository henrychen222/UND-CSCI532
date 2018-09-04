#!/bin/sh
mpirun -np 5 kmeans_mpi --cluster_num 3 --star_files ./stars/stars-9.txt ./stars/stars-15.txt ./stars/stars-82.txt
