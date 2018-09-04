1. compile 
   nvcc kernel.cu -o ass4

2. use "nvprof" command to run by output running time
   nvprof ./ass4

3. check the running result and output time printed by nvprof  (Has Avg,Max and Min time for "convKernel()", "CUDA memcpy DtoH" and "CUDA memcpy HtoD" respectively)

4. change "srcHeight", "srcWidth" and "MaskHeight", "MaskWeight" to 
       32x32, 64x64, 128x128, 256x256, 512x512, 1024x1024 (input)
       4x4, 8x8, 16x16, 32x32, 64x64 and 128x128  (filter)
    and calculate the average runtime of each combination in a Table(check Excel or PDF file)
