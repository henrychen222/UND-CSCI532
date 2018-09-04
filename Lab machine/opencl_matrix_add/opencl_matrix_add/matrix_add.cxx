#include <cmath>
using std::sqrt;

#include <chrono>

#include <iomanip>
using std::fixed;
using std::setw;
using std::setprecision;

#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

#ifdef __OPENCL__

#include <cstdio>
#include <cstring>
#include <time.h>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "opencl_utils.hxx"

#endif


//set up global opencl variables
cl_device_id device;
cl_context context;
cl_program matrix_add_kernel;
cl_kernel kernel;
cl_command_queue queue;

cl_mem sizes_oepncl;
cl_mem A_opencl;
cl_mem B_opencl;
cl_mem C_opencl;

size_t *global_size;
size_t *local_size;

#define MATRIX_ADD_KERNEL_FILE "../matrix_add_kernel.cl"


//initialize the memory for a 3D array
void initialize_3d(float ****v, uint32_t v_z, uint32_t v_y, uint32_t v_x) {
    (*v) = (float***)malloc(sizeof(float**) * v_z);
    for (uint32_t z = 0; z < v_z; z++) {
        (*v)[z] = (float**)malloc(sizeof(float*) * v_y);
        for (uint32_t y = 0; y < v_y; y++) {
            (*v)[z][y] = (float*)malloc(sizeof(float) * v_x);
        }
    }
}

//set the values in a 3d array to random numbers
void set_to_random_3d(float ***v, uint32_t v_z, uint32_t v_y, uint32_t v_x) {
    for (uint32_t z = 0; z < v_z; z++) {
        for (uint32_t y = 0; y < v_y; y++) {
            for (uint32_t x = 0; x < v_x; x++) {
                v[z][y][x] = drand48();
            }
        }
    }
}

//set the values in a 3d array to 0
void set_to_zero_3d(float ***v, uint32_t v_z, uint32_t v_y, uint32_t v_x) {
    for (uint32_t z = 0; z < v_z; z++) {
        for (uint32_t y = 0; y < v_y; y++) {
            for (uint32_t x = 0; x < v_x; x++) {
                v[z][y][x] = 0.0;
            }
        }
    }
}

//copy the values from a 3d array to a flattened 1d array
void copy_3d_to_1d(float ***input, uint32_t input_z, uint32_t input_y, uint32_t input_x, float *output) {
    uint32_t current_output = 0;
    for (uint32_t z = 0; z < input_z; z++) {
        for (uint32_t y = 0; y < input_y; y++) {
            for (uint32_t x = 0; x < input_x; x++) {
                output[current_output++] = input[z][y][x];
            }
        }
    }
}

//copy the values from a flattened 1d array to a 3d array
void copy_1d_to_3d(float *input, uint32_t output_z, uint32_t output_y, uint32_t output_x, float ***output) {
    uint32_t current_output = 0;
    for (uint32_t z = 0; z < output_z; z++) {
        for (uint32_t y = 0; y < output_y; y++) {
            for (uint32_t x = 0; x < output_x; x++) {
                output[z][y][x] = input[current_output++];
            }
        }
    }
}

void print_3d(string name, float ***input, uint32_t input_z, uint32_t input_y, uint32_t input_x) {
    cout << "MATRIX '" << name << "'" << endl;
    for (uint32_t z = 0; z < input_z; z++) {
        for (uint32_t y = 0; y < input_y; y++) {
            for (uint32_t x = 0; x < input_x; x++) {
                cout << setw(10) << fixed << input[z][y][x];
            }
            cout << endl;
        }
        cout << endl;
    }
}


void matrix_add(float ***A, float ***B, float ***C, int input_z, int input_y, int input_x) {
    for (uint32_t z = 0; z < input_z; z++) {
        for (uint32_t y = 0; y < input_y; y++) {
            for (uint32_t x = 0; x < input_x; x++) {
                C[z][y][x] = A[z][y][x] + B[z][y][x];
            }
        }
    }
 }

bool equal_3d(float ***M1, float ***M2, int input_z, int input_y, int input_x) {
    for (uint32_t z = 0; z < input_z; z++) {
        for (uint32_t y = 0; y < input_y; y++) {
            for (uint32_t x = 0; x < input_x; x++) {
                if (M1[z][y][x] != M2[z][y][x]) return false;
            }
        }
    }
    return true;
}

void initialize_opencl(int z, int y, int x) {
    //OpenCL structures
    cl_int err;

    //Create device and context
    device = create_device();
    context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    check_error(err, "couldn't create a context, err: %d", err);

    //Create a command queue
    queue = clCreateCommandQueue(context, device, 0, &err);
    check_error(err, "couldn't create a command queue: %d", err);

    // Build program
    matrix_add_kernel = build_program(context, device, MATRIX_ADD_KERNEL_FILE);

    // Create a kernel
    kernel = clCreateKernel(matrix_add_kernel, "matrix_add", &err);
    check_error(err, "couldn't create a kernel: %d", err);


    //A_opencl, B_opencl, and C_opencl are set as global variables so we can reuse them
    int size = sizeof(float) * z * y * x;
    A_opencl = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &err);
    check_error(err, "could not create A_opencl buffer: %d", err);

    B_opencl = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &err);
    check_error(err, "could not create B_opencl buffer: %d", err);

    C_opencl = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &err);
    check_error(err, "could not create C_opencl buffer: %d", err);

    // only need to set the kernel arguments once, and we can then re-use them if
    // those cl_mem variables don't change
    // Create kernel arguments
    err =  clSetKernelArg(kernel, 0, sizeof(cl_mem), &A_opencl);
    check_error(err, "couldn't create A_opencl argument: %d", err);

    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &B_opencl);
    check_error(err, "couldn't create B_opencl argument: %d", err);

    err = clSetKernelArg(kernel, 2, sizeof(cl_mem), &C_opencl);
    check_error(err, "couldn't create C_opencl argument: %d", err);

    global_size = (size_t*)malloc(sizeof(size_t) * 3);
    global_size[0] = z;
    global_size[1] = y;
    global_size[2] = x;

    local_size = (size_t*)malloc(sizeof(size_t) * 3);
    local_size[0] = 5;
    local_size[1] = 5;
    local_size[2] = 5;
}


void matrix_add_opencl(float *A_flat, float *B_flat, float *C_flat, int z, int y, int x) {
    cl_int err;

    int size = sizeof(float) * z * y * x;

    err = clEnqueueWriteBuffer(queue, A_opencl, CL_TRUE, 0, size, A_flat, 0, NULL, NULL);
    check_error(err, "couldn't write to the A_opencl buffer: %d", err);

    err = clEnqueueWriteBuffer(queue, B_opencl, CL_TRUE, 0, size, B_flat, 0, NULL, NULL);
    check_error(err, "couldn't write to the A_opencl buffer: %d", err);

    // Enqueue kernel
    err = clEnqueueNDRangeKernel(queue, kernel, 3, NULL, global_size, local_size, 0, NULL, NULL); 
    check_error(err, "couldn't enqueue the kernel: %d", err);

    err = clFinish(queue);
    check_error(err, "queue errored on finish: %d", err);

    // Read the kernel's output
    err = clEnqueueReadBuffer(queue, C_opencl, CL_TRUE, 0, size, C_flat, 0, NULL, NULL);
    check_error(err, "couldn't read the C_opencl buffer: %d", err);
}



int main(int argc, char **argv) {
    vector<string> arguments = vector<string>(argv, argv + argc);

    uint32_t x = 10;
    uint32_t y = 10;
    uint32_t z = 10;

    float ***A;
    float ***B;
    float ***C;
    float ***C_gpu;

    initialize_3d(&A, z, y, x);
    initialize_3d(&B, z, y, x);
    initialize_3d(&C, z, y, x);
    initialize_3d(&C_gpu, z, y, x);

    set_to_random_3d(A, z, y, x);
    set_to_random_3d(B, z, y, x);
    set_to_zero_3d(C, z, y, x);
    set_to_zero_3d(C_gpu, z, y, x);

    //print_3d("A", A, z, y, x);
    //print_3d("B", B, z, y, x);

    float *A_flat = (float*)malloc(sizeof(float) * z * y * x);
    float *B_flat = (float*)malloc(sizeof(float) * z * y * x);
    float *C_flat = (float*)malloc(sizeof(float) * z * y * x);

    copy_3d_to_1d(A, z, y, x, A_flat);
    copy_3d_to_1d(B, z, y, x, B_flat);
    copy_3d_to_1d(C, z, y, x, C_flat);
    cout << "created initial arrays." << endl;

    cout << endl << "initializing opencl." << endl;
    initialize_opencl(z, y, x);
    cout << "initialized successfully." << endl;

    using namespace std::chrono;

    high_resolution_clock::time_point t1, t2;
    duration<float, std::milli> time_span;

    t1 = high_resolution_clock::now();
    matrix_add_opencl(A_flat, B_flat, C_flat, z, y, x);
    t2 = high_resolution_clock::now();

    time_span = t2 - t1;

    cout << "OpenCL Matrix Add took: " << time_span.count() / 1000.0 << " seconds." << endl << endl;

    copy_1d_to_3d(C_flat, z, y, x, C_gpu);
    //print_3d("C_GPU", C_gpu, z, y, x);


    t1 = high_resolution_clock::now();
    matrix_add(A, B, C, z, y, x);
    t2 = high_resolution_clock::now();

    time_span = t2 - t1;

    cout << "CPU Matrix Add took: " << time_span.count() / 1000.0 << " seconds." << endl << endl;

    //print_3d("C_CPU", C, z, y, x);

    cout << "Matrices equal? " << equal_3d(C_gpu, C, z, y, x) << endl;

    clReleaseMemObject(A_opencl);
    clReleaseMemObject(B_opencl);
    clReleaseMemObject(C_opencl);

    clReleaseKernel(kernel);
    clReleaseCommandQueue(queue);
}


