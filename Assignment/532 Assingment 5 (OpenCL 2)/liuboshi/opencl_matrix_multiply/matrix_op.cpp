#include <cmath>
using std::sqrt;

#include <chrono>

#include <iomanip>
using std::fixed;
using std::setw;
using std::setprecision;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#define MATRIX_OP_KERNEL_FILE "../matrix_multiply.cl"

void check_error(cl_int err, const char* fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);

    if (err < 0)
    {
        vfprintf(stderr, fmt, argp);
        fprintf(stderr, "\n");
        exit(1);
    };
}

//initialize the memory for a 3D array
void initialize_2d(float ***v, int row, int col)
{
    (*v) = new float*[row];
    for (int idx = 0; idx < row; ++idx)
    {
        (*v)[idx] = new float[col];
    }
}

//set the values in a 3d array to random numbers
void set_to_random_2d(float **v, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v[i][j] = static_cast<float>(drand48 ());
        }
    }
}

//set the values in a 3d array to 0
void set_to_zero_2d(float **v, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v[i][j] = 0.0f;
        }
    }
}

//copy the values from a 3d array to a flattened 1d array
void copy_2d_to_1d(float ** matrix2d, int row, int col, float * array)
{
    int pos = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            array[pos++] = matrix2d[i][j];
        }
    }
}

//copy the values from a flattened 1d array to a 3d array
void copy_1d_to_2d(float * array, int row, int col, float ** matrix2d)
{
    int pos = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            matrix2d[i][j] = array[pos++];
        }
    }
}

void print_2d(string name, float ** matrix2d, int row, int col)
{
    cout << "MATRIX '" << name << "'" << endl;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            cout << setw(8) << setprecision (3) << fixed << matrix2d[i][j];
        }
        cout << endl;
    }
}

bool equal_2d(float **M1, float **M2, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (M1[i][j] - M2[i][j] > 0.0001f || M1[i][j] - M2[i][j] < -0.0001f)
            {
                //cerr << setprecision (10) << "M1 = " << M1[i][j] << "  M2 = " << M2[i][j] << endl;
                return false;
            }
        }
    }
    return true;
}

int main(int argc, char **argv)
{
    //step 0: prepare program
    int m = 0, n = 0, p = 0, q = 0;
    int argset = 0;

    //parse command line parameters
    for (int idx = 0; idx < argc; ++idx)
    {
        if (0 == strcmp(argv[idx], "-m"))
        {
            m = atoi(argv[++idx]);
            argset++;
            continue;
        }
        if (0 == strcmp(argv[idx], "-n"))
        {
            n = atoi(argv[++idx]);
            argset++;
            continue;
        }
        if (0 == strcmp(argv[idx], "-p"))
        {
            p = atoi(argv[++idx]);
            argset++;
            continue;
        }
        if (0 == strcmp(argv[idx], "-q"))
        {
            q = atoi(argv[++idx]);
            argset++;
            continue;
        }
    }

    if (4 != argset)
    {
        cerr << "No Proper arguments!" << endl;
        exit(1);
    }

    if (n != p)
    {
        cerr << "The multiply operation is undefined between two matrix: Matrix(" << m << "," << n << ") multiply Matrix(" << p << "," << q << ")." << endl;
        exit(1);
    }

    //Initialize matrixs
    float **A;
    float **B;
    float **C;
    float **C_gpu;

    initialize_2d(&A, m, n);
    initialize_2d(&B, p, q);
    initialize_2d(&C, m, q);
    initialize_2d(&C_gpu, m, q);

    set_to_random_2d(A, m, n);
    set_to_random_2d(B, p, q);
    set_to_zero_2d(C, m, q);
    set_to_zero_2d(C_gpu, m, q);

    float * A_flat = new float[m*n];
    float * B_flat = new float[p*q];
    float * C_flat = new float[m*q];

    copy_2d_to_1d(A, m, n, A_flat);
    copy_2d_to_1d(B, p, q, B_flat);

    print_2d("A", A, m, n);
    print_2d("B", B, p, q);

    //Step 1: run a CPU 2D matrix multiply
    /// NOTE Question 1: 2d multiply
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<float, std::milli> time_span;
    t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < q; j++)
        {
            C[i][j] = 0.0f;
            for (int rang = 0; rang < n; rang++)
            {
                C[i][j] += A[i][rang]*B[rang][j];
            }
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span = t2 - t1;
    cout << "CPU Matrix Mutiply took: " << time_span.count() / 1000.0f << " seconds." << endl << endl;
    /// End of Question1

    //Step 2: prepare for an OpenCL run
    cout << endl << "initializing opencl." << endl;
    cl_int err;
    cl_device_id device;
    cl_context context;
    cl_program matrix_program;
    cl_kernel kernel1, kernel2;
    cl_command_queue queue;
    cl_mem A_opencl;
    cl_mem B_opencl;
    cl_mem C_opencl;
    int C_Row = m;
    int C_Col = q;
    int C_Num = n;

    size_t global_size[2] = {static_cast<size_t>(C_Row), static_cast<size_t>(C_Col)};
    //WARING: Adjust this to get over some opencl runtime error correspond to group_size
    size_t local_size[2] = {16, 16};

    //create device and context
    cl_platform_id platform;
    cl_device_id dev[3];
    cl_uint num_devices;

    //identify a platform
    err = clGetPlatformIDs(1, &platform, NULL);
    check_error(err, "Couldn't identify a platform");

    //access a device
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 3, dev, &num_devices);
    if (err == CL_DEVICE_NOT_FOUND)
    {
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 3, dev, &num_devices);
    }

    check_error(err, "Couldn't access any devices");

    char deviceName[1024];              //this string will hold the devices name
    char vendor[1024];                  //this strirng will hold a platforms vendor
    cl_uint numberOfCores;              //this variable holds the number of cores of on a device
    cl_long amountOfMemory;             //this variable holds the amount of memory on a device
    cl_uint clockFreq;                  //this variable holds the clock frequency of a device
    cl_ulong maxAlocatableMem;          //this variable holds the maximum allocatable memory
    cl_ulong localMem;                  //this variable holds local memory for a device
    cl_bool available;                  //this variable holds if the device is available
    cl_ulong constantBufferSize;
    cl_uint maxWorkItemDimensions;
    size_t maxWorkItemSizes[3];
    size_t maxWorkGroupSize;

    for (uint32_t i = 0; i < num_devices; i++) {
        //scan in device information
        clGetDeviceInfo(dev[i], CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_VENDOR, sizeof(vendor), vendor, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(numberOfCores), &numberOfCores, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(amountOfMemory), &amountOfMemory, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clockFreq), &clockFreq, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(maxAlocatableMem), &maxAlocatableMem, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localMem), &localMem, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_AVAILABLE, sizeof(available), &available, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(constantBufferSize), &constantBufferSize, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(maxWorkGroupSize), &maxWorkGroupSize, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(maxWorkItemDimensions), &maxWorkItemDimensions, NULL);
        clGetDeviceInfo(dev[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(maxWorkItemSizes), &maxWorkItemSizes, NULL);

        //print out device information
        printf("\tDevice:\n");
        printf("\t\tName:\t\t\t\t%s\n", deviceName);
        printf("\t\tVendor:\t\t\t\t%s\n", vendor);
        printf("\t\tAvailable:\t\t\t%s\n", available ? "Yes" : "No");
        printf("\t\tCompute Units:\t\t\t%u\n", numberOfCores);
        printf("\t\tClock Frequency:\t\t%u mHz\n", clockFreq);
        printf("\t\tGlobal Memory:\t\t\t%0.00f mb\n", (double)amountOfMemory/1048576);
        printf("\t\tMax Allocateable Memory:\t%0.00f mb\n", (double)maxAlocatableMem/1048576);
        printf("\t\tLocal Memory:\t\t\t%u kb\n\n", (unsigned int)localMem);
        printf("\t\tMax Constant Buffer Size:\t%0.00f mb\n\n", (double)constantBufferSize/1048576);
        printf("\t\tMax Work Group Size:\t\t%zu\n", maxWorkGroupSize);
        printf("\t\tMax Work Item Dimensions:\t%u\n", maxWorkItemDimensions);
        printf("\t\tMax Work Item Sizes:\t\t%zu %zu %zu\n", maxWorkItemSizes[0], maxWorkItemSizes[1], maxWorkItemSizes[2]);
        printf("\n");
    }

    device = dev[0];

    context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    check_error(err, "couldn't create a context, err: %d", err);

    //Create a command queue
    queue = clCreateCommandQueue(context, device, 0, &err);
    check_error(err, "couldn't create a command queue: %d", err);

    //Read program file and place content into buffer
    FILE * program_handle = fopen(MATRIX_OP_KERNEL_FILE, "r");
    if (program_handle == NULL)
    {
        fprintf(stderr, "Couldn't find the program file: '%s'\n", MATRIX_OP_KERNEL_FILE);
        return 1;
    }
    fseek(program_handle, 0, SEEK_END);
    size_t program_size = (size_t)ftell(program_handle);
    rewind(program_handle);
    char * program_buffer = new char[program_size + 1];
    program_buffer[program_size] = '\0';
    fread(program_buffer, sizeof(char), program_size, program_handle);
    fclose(program_handle);

    //Create program from file
    matrix_program = clCreateProgramWithSource(context, 1, (const char**)&program_buffer, &program_size, &err);
    free(program_buffer);

    //build program
    err = clBuildProgram(matrix_program, 0, NULL, NULL, NULL, NULL);
    check_error(err, "Couldn't build program: '%d'", err);

    //create a kernels
    kernel1 = clCreateKernel(matrix_program, "matrix_multiply", &err);
    check_error(err, "couldn't create a kernel: %d", err);

    //create buffer and copy data. use global data so that the matrix can be larger in constant memory!
    A_opencl = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, m*n*sizeof(float), A_flat, &err);
    check_error(err, "could not create A_opencl buffer: %d", err);

    B_opencl = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, p*q*sizeof(float), B_flat, &err);
    check_error(err, "could not create B_opencl buffer: %d", err);

    C_opencl = clCreateBuffer(context, CL_MEM_READ_WRITE, m*q*sizeof(float), NULL, &err);
    check_error(err, "could not create C_opencl buffer: %d", err);

    //Step 3: config and run kernel : matrix_multiply with no tiling
    ///NOTE: Question 2 OpenCL matrix multiply with no tiling
    //config kernel1 and start the kernel.
    err = clSetKernelArg(kernel1, 0, sizeof(cl_mem), &A_opencl);
    check_error(err, "couldn't create A_opencl argument: %d", err);

    err = clSetKernelArg(kernel1, 1, sizeof(cl_mem), &B_opencl);
    check_error(err, "couldn't create B_opencl argument: %d", err);

    err = clSetKernelArg(kernel1, 2, sizeof(cl_mem), &C_opencl);
    check_error(err, "couldn't create C_opencl argument: %d", err);

    err = clSetKernelArg(kernel1, 3, sizeof(int), &C_Num);
    check_error(err, "couldn't create C_Row argument: %d", err);

    cout << "kernel(matrix_multiply(...)) initialized successfully." << endl;

    t1 = std::chrono::high_resolution_clock::now();
    //enqueue kernel
    err = clEnqueueNDRangeKernel(queue, kernel1, 2, NULL, global_size, local_size, 0, NULL, NULL);
    check_error(err, "couldn't enqueue the kernel: %d", err);
    err = clFinish(queue);
    check_error(err, "queue errored on finish: %d", err);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = t2 - t1;
    cout << endl << "OpenCL Matrix Multiply by kernel(matrix_multiply(...)) took: " << time_span.count() / 1000.0f << " seconds." << endl << endl;

    //read the kernel's output
    err = clEnqueueReadBuffer(queue, C_opencl, CL_TRUE, 0, static_cast<size_t>(m*q*sizeof(float)), C_flat, 0, NULL, NULL);
    check_error(err, "couldn't read the C_opencl buffer: %d", err);

    copy_1d_to_2d(C_flat, m, q, C_gpu);

    //print_2d("C_CPU", C, m, q);
    //print_2d("C_GPU", C_gpu, m, q);

    cout << "Matrices equal? " << (equal_2d(C_gpu, C, m, q) ? "Yes" : "No") << endl;
    /// End of Question2

    //Step 4: config and run kernel with tiling
    /// NOTE: Question3 OpenCL matrix multiply with tiling 16
    kernel2 = clCreateKernel(matrix_program, "matrix_multiply_tiling", &err);
    check_error(err, "Couldn't create a kernel: %d", err);

    err = clSetKernelArg(kernel2, 0, sizeof(cl_mem), &A_opencl);
    check_error(err, "couldn't create A_opencl argument: %d", err);

    err = clSetKernelArg(kernel2, 1, sizeof(cl_mem), &B_opencl);
    check_error(err, "couldn't create B_opencl argument: %d", err);

    err = clSetKernelArg(kernel2, 2, sizeof(cl_mem), &C_opencl);
    check_error(err, "couldn't create C_opencl argument: %d", err);

    err = clSetKernelArg(kernel2, 3, sizeof(int), &C_Num);
    check_error(err, "couldn't create C_Row argument: %d", err);

    cout << "kernel(matrix_multiply_tiling(...)) initialized successfully." << endl;

    t1 = std::chrono::high_resolution_clock::now();
    // Enqueue kernel
    err = clEnqueueNDRangeKernel(queue, kernel2, 2, NULL, global_size, local_size, 0, NULL, NULL);
    check_error(err, "couldn't enqueue the kernel: %d", err);
    err = clFinish(queue);
    check_error(err, "queue errored on finish: %d", err);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = t2 - t1;
    cout << endl << "OpenCL Matrix Multiply by kernel(matrix_multiply(...)) took: " << time_span.count() / 1000.0f << " seconds." << endl << endl;

    // Read the kernel's output
    err = clEnqueueReadBuffer(queue, C_opencl, CL_TRUE, 0, static_cast<size_t>(m*q*sizeof(float)), C_flat, 0, NULL, NULL);
    check_error(err, "couldn't read the C_opencl buffer: %d", err);

    copy_1d_to_2d(C_flat, m, q, C_gpu);

    //print_2d("C_CPU", C, m, q);
    //print_2d("C_GPU", C_gpu, m, q);

    cout << "Matrices equal between CPU and NO-TILING-KERNEL? " << (equal_2d(C_gpu, C, m, q) ? "Yes" : "No") << endl;
    /// End of Question3

    delete [] A_flat;
    delete [] B_flat;
    delete [] C_flat;

    for (int i = 0; i < m; ++i)
    {
        delete [] A[i];
    }
    delete [] A;

    for (int i = 0; i < p; ++i)
    {
        delete [] B[i];
    }
    delete [] B;

    for (int i = 0; i < m; ++i)
    {
        delete [] C[i];
        delete [] C_gpu[i];
    }
    delete [] C;
    delete [] C_gpu;

    clReleaseMemObject(A_opencl);
    clReleaseMemObject(B_opencl);
    clReleaseMemObject(C_opencl);

    clReleaseKernel(kernel1);
    clReleaseProgram (matrix_program);
    clReleaseCommandQueue(queue);
    clReleaseContext (context);

    return 0;
}


