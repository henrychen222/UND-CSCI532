
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const int srcHeight = 32;
const int srcWidth = 32;
const int maskHeight = 4;
const int maskWidth = 4;
const int dstHeight = srcHeight - maskHeight + 1;
const int dstWidth = srcWidth - maskWidth + 1;

const int BLOCK_DIM_X = 32;
const int BLOCK_DIM_Y = 32;

__constant__ float devMask[maskHeight*maskWidth];

__global__ void convKernel(int maskW, int maskH, int srcW, int srcH, int dstW, int dstH, float * src, float * dst)
{
    const int tx = blockDim.x * blockIdx.x + threadIdx.x;
	const int ty = blockDim.y * blockIdx.y + threadIdx.y;

	//coordination outside the range, drop them
	if (tx < dstW && ty < dstH)
	{
		dst[ty*dstW + tx] = 0.0;
		for (int cy = 0; cy < maskH; cy++)
		{
			for (int cx = 0; cx < maskW; cx++)
			{
				dst[ty*dstW + tx] += devMask[cy*maskW + cx] * src[ty*srcW + tx];
			}
		}
	}

	return;
}

void printMatrix(const char * name, float * matrix, int height, int width)
{
	printf("%s\n", name);
	for (int idx = 0; idx < height*width; ++idx)
	{
		printf("%4.1f ", matrix[idx]);
		if (idx % width == width - 1)
		{
			printf("\n");
		}
	}
	printf("\n");
}

int main()
{
	//list and select device
	cudaDeviceProp prop;
	int devCount = 0;
	cudaError_t cudaError = cudaGetDeviceCount(&devCount);
	for (int i = 0; i < devCount; ++i)
	{
		cudaError = cudaGetDeviceProperties(&prop, i);
		printf("Cuda Device %d/%d:\n", i, devCount);
		printf("Cuda Device Name: %s\n", prop.name);
		printf("Global Memory: %.2lf MB\n", prop.totalGlobalMem / 1024.0 / 1024.0);
		printf("Shared Memory Per Block: %.2lf KB\n", prop.sharedMemPerBlock / 1024.0);
		printf("Register Per Block: %d\n", prop.regsPerBlock);
		printf("Max Thread Per Block: %d\n", prop.maxThreadsPerBlock);
		printf("Max Size of each dim of block xyz[%d %d %d]\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
		printf("Compute Capability: %d.%d\n", prop.major, prop.minor);
		
	}
	int selectedDev = 0;
	while (1)
	{
		printf("Select a device by device ID: ");
		scanf("%d", &selectedDev);
		if (cudaSuccess != cudaSetDevice(selectedDev))
		{
			printf("\nInvalid input, try again!\n");
		}
		else
		{
			printf("Device #%d is now in use!\n", selectedDev);
			break;
		}
	}

	//allocate host memory
	float * srcMatrix = (float *)calloc(srcHeight*srcWidth, sizeof(float));
	float * dstMatrix = (float *)calloc(dstHeight*dstWidth, sizeof(float));
	float * maskMatrix = (float *)calloc(maskHeight*maskWidth, sizeof(float));

	for (int idx = 0; idx < srcHeight*srcWidth; ++idx)
	{
		//srand((unsigned)time(NULL) + idx);
		srcMatrix[idx] = 1.1;//(rand() % 100) / 100.0;
	}

	printMatrix("srcMatrix", srcMatrix, srcHeight, srcWidth);

	for (int idx = 0; idx < maskHeight*maskWidth; ++idx)
	{
		//srand((unsigned)time(NULL) + idx);
		maskMatrix[idx] = 1.0;//(rand() % 100) / 100.0;
	}

    //allocate device memory
	float * devSrc;
	float * devDst;

	cudaMalloc((float **)&devSrc, srcHeight*srcWidth * sizeof(float));
	cudaMalloc((float **)&devDst, dstHeight*dstWidth * sizeof(float));

	cudaMemcpyToSymbol(devMask, maskMatrix, maskHeight*maskWidth * sizeof(float));
	cudaMemcpy(devSrc, srcMatrix, srcHeight*srcWidth * sizeof(float), cudaMemcpyHostToDevice);

	//call kernel to do parallel compute
	dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
	dim3 grid((srcWidth + block.x - 1) / block.x, (srcHeight + block.y - 1) / block.y);
	convKernel << <grid, block >> > (maskWidth, maskHeight, srcWidth, srcHeight, dstWidth, dstHeight, devSrc, devDst);
	cudaDeviceSynchronize();

	cudaMemcpy(dstMatrix, devDst, dstHeight*dstWidth * sizeof(float), cudaMemcpyDeviceToHost);

	printMatrix("dstMatrix", dstMatrix, dstHeight, dstWidth);

	cudaFree(devDst);
	cudaFree(devSrc);

	free(srcMatrix);
	free(dstMatrix);
	free(maskMatrix);

    return 0;
}
