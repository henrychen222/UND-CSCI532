#!/bin/sh
mpirun -np 5 ass2 --chr_files ./chr1.fa ./chr2.fa --seed_files ./ERR192339.fastq ./ERR192338.fastq
