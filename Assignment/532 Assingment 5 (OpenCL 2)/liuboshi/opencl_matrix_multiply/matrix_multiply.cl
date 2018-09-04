#pragma unroll

__kernel void matrix_multiply(__global float * A, __global float * B, __global float * C, const int P)
{
    int row = get_global_id(0);
    int col = get_global_id(1);

    int COL = get_global_size(1);

    float tmp = 0.0;
    for (int k = 0; k < P; ++k)
    {
        tmp += A[row * P + k] * B[k * COL + col];
    }
    C[row * COL + col] = tmp;
}

__kernel void matrix_multiply_tiling(__global float * A, __global float * B, __global float * C, const int P)
{
    int local_row_id = get_local_id(0);
    int local_col_id = get_local_id(1);

    int global_row_id = get_global_id(0);
    int global_col_id = get_global_id(1);

    int group_row_id = get_group_id(0);
    int group_col_id = get_group_id(1);

    __local float Ashare[16][16];
    __local float Bshare[16][16];

    float result = 0.0f;

    for (int p = 0; p < P / 16; ++p)
    {
        //load all shared mem
        Ashare[local_row_id][local_col_id] = A[global_row_id*P+p*16+local_col_id];
        Bshare[local_row_id][local_col_id] = B[(p*16+local_row_id)*P+global_col_id];
        barrier(CLK_LOCAL_MEM_FENCE);

        for (int i = 0; i < 16; ++i)
        {
            result += Ashare[local_row_id][i]*Bshare[i][local_col_id];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    C[global_row_id * P + global_col_id] = result;
}
