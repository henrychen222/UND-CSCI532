1. mkdir - p build
2. cd build
3. cmake ..
4. make
5. ./matrix_op -m 512 -n 512 -p 512 -q 512 (m n p q should be interger times of 16 using opencl1.x or opecl runtime will report an error)
