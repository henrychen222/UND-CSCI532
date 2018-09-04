#include <string.h>
#include <iomanip>
using std::fixed;
using std::setprecision;
using std::setw;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ofstream;

#include <sstream>
using std::ostringstream;

#include <string>
using std::string;

#include <mpi/mpi.h>

//write the state of the smulation to a file
void write_simulation_state(string name, int height, int width, int time_step, double **values) {
    ostringstream filename;
    filename << name << "_" << height << "x" << width << "_" << time_step;

    ofstream outfile(filename.str());

    for (uint32_t i = 0; i < height; i++) {
        for (uint32_t j = 0; j < width; j++) {
            outfile << setw(10) << fixed << setprecision(5) << values[i][j];
        }
        outfile << endl;
    }
}

int main(int argc, char **argv) {

    //mpi init
    MPI_Init(&argc, &argv);
    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    string simulation_name;
    //the height and width of the simulation
    int height;
    int width;
    //horizontal and vertical slices will be used to determine how to divide up your MPI processes
    int vertical_slices;
    int horizontal_slices;
    //how long to run the simulation for
    int time_steps;

    for(int idx = 1; idx < argc; ++idx) {
        if (0 == strcmp(argv[idx], "--sim_name")) simulation_name.append(argv[++idx]);
        if (0 == strcmp(argv[idx], "--height")) height = atoi(argv[++idx]);
        if (0 == strcmp(argv[idx], "--width")) width = atoi(argv[++idx]);
        if (0 == strcmp(argv[idx], "--v_slice")) vertical_slices = atoi(argv[++idx]);
        if (0 == strcmp(argv[idx], "--h_slice")) horizontal_slices = atoi(argv[++idx]);
        if (0 == strcmp(argv[idx], "--time_step")) time_steps = atoi(argv[++idx]);
    }

    cout << "Proc " << comm_rank << " :"\
         << " sim_name=" << simulation_name
         << " height=" << height
         << " width=" << width
         << " v_slice=" << vertical_slices
         << " h_slice=" << horizontal_slices
         << " time_step=" << time_steps << endl;

    //You will need to have a number of MPI processes equal to vertical_slices * horizontal_slices
    double *values = NULL;
    double **values_reo = NULL;

    MPI_Status status[4];

    if (0 == comm_rank) {
        values = new double[height*width];
        for (int i = 0; i < height*width; ++i) values[i] = 0.0;
        values_reo = new double*[height + 2];
        for (int i = 0; i < height + 2; ++i) {
            values_reo[i] = new double[width + 2];
            for (int j = 0; j < width + 2; ++j) values_reo[i][j] = 0.0;
        }
        //put heat source
        for (int i = 0; i < height + 2; ++i) {
            values_reo[i][0] = 1.0;
            values_reo[i][width+1] = -1.0;
        }
    }

    //slice locate
    const int row_num = comm_rank / vertical_slices;//row of this rectangle block
    const int col_num = comm_rank % vertical_slices;//col of this rectangle block
    const int row_max = vertical_slices;
    const int col_max = horizontal_slices;

    //compute neighbor ranks, invalid neighbor rank will be rank of itself
    int west_rank = row_num * row_max + (col_num - 1);
    if (0 == col_num) west_rank = -1;
    int north_rank = (row_num - 1) * row_max + col_num;
    if (0 == row_num) north_rank = -1;
    int east_rank = row_num * row_max + (col_num + 1);
    if (col_max - 1 == col_num) east_rank = -1;
    int south_rank = (row_num + 1) * row_max + col_num;
    if (row_max - 1 == row_num) south_rank = -1;

    //allocate local data and init
    const uint32_t my_height = height / vertical_slices;
    const uint32_t my_width = width / horizontal_slices;
    double ** my_block = new double*[my_height + 2];
    for (uint32_t i = 0; i < my_height + 2; ++i) {
        my_block[i] = new double[my_width + 2];
        for (uint32_t j = 0; j < my_width + 2; ++j) {
            my_block[i][j] = 0.0;
        }
    }

    //allocate result array
    double my_result[my_height*my_width];
    for (uint32_t i = 0; i < my_height*my_width; ++i) my_result[i] = 0.0;

    //allocate halo array for transfer
    double from_east[my_height];
    double to_east[my_height];
    double from_west[my_height];
    double to_west[my_height];
    double from_north[my_width];
    double to_north[my_width];
    double from_south[my_width];
    double to_south[my_width];

    //put heat source
    for (int i = 0; i < my_height; ++i) from_west[i] = 0.0;
    for (int i = 0; i < my_height; ++i) from_east[i] = 0.0;
    for (int i = 0; i < my_width; ++i) from_south[i] = 0.0;
    for (int i = 0; i < my_width; ++i) from_north[i] = 0.0;
    if (0 == col_num) {
        for (int i = 0; i < my_height; ++i) {
            from_west[i] = 1.0;
        }
    }
    if (col_max - 1 == col_num) {
        for (int i = 0; i < my_height; ++i) {
            from_east[i] = -1.0;
        }
    }

    //after intial halo transfer, start to iteration
    for (int steps = 0; steps < time_steps; steps++)
    {
        //fill in halo with updates
        for (int i = 0; i < my_height; ++i) my_block[i + 1][0] = from_west[i];
        for (int i = 0; i < my_height; ++i) my_block[i + 1][my_width + 1] = from_east[i];
        for (int i = 0; i < my_width; ++i) my_block[my_height + 1][i + 1] = from_south[i];
        for (int i = 0; i < my_width; ++i) my_block[0][i + 1] = from_north[i];

        //perform a calculation of the heat diffusion
        for (uint32_t i = 0; i < my_height; ++i) {
            for (uint32_t j = 0; j < my_width; ++j) {
                my_result[i*my_width+j] = \
                        (my_block[i+1][j] + my_block[i+1][j+2] + my_block[i][j+1] + my_block[i+2][j+1])/4.0;
            }
        }

        //update local data
        for (int i = 0; i < my_height; ++i) {
            for (int j = 0; j < my_width; ++j) {
                my_block[i + 1][j + 1] = my_result[i * my_width + j];
            }
        }

        /*
        //WARNING: write debug files, do not open
        char arr[128] = {'\0'};
        char arr1[2] = {48+comm_rank, '\0'};
        strcat(arr, simulation_name.c_str());
        strcat(arr, arr1);
        string name(arr);
        write_simulation_state(name, my_height+2, my_width+2, steps, my_block);
        */

        //Gather blocks to rank 0, then reo
        MPI_Gather(my_result, my_width*my_height, MPI_DOUBLE, values, my_width*my_height, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (0 == comm_rank) {
            for (int rk = 0; rk < col_max*row_max; ++rk) {
                int block_row = rk / col_max;
                int block_col = rk % col_max;
                for (int idx_row = 0; idx_row < my_height; ++idx_row) {
                    for (int idx_col = 0; idx_col < my_width; ++idx_col) {
                        values_reo[block_row * my_height + idx_row + 1][block_col * my_width + idx_col + 1] = \
                                values[rk * my_height * my_width + idx_row*my_width + idx_col];
                    }
                }
            }

            write_simulation_state(simulation_name, height+2, width+2, steps, values_reo);
        }

        //extract halo to transfer
        for (int i = 0; i < my_height; ++i) to_east[i] = my_result[i * my_width + (my_width - 1)];
        for (int i = 0; i < my_height; ++i) to_west[i] = my_result[i * my_width + 0];
        for (int i = 0; i < my_width; ++i) to_north[i] = my_result[0 * my_width + i];
        for (int i = 0; i < my_width; ++i) to_south[i] = my_result[(my_height - 1) * my_width  + i];

        //halo tansfer via MPI
        //east<->west boarder exchange
        if (east_rank != -1) MPI_Send(to_east, my_height, MPI_DOUBLE, east_rank, 1, MPI_COMM_WORLD);
        if (west_rank != -1) MPI_Recv(from_west, my_height, MPI_DOUBLE, west_rank, 1, MPI_COMM_WORLD, &status[0]);
        if (west_rank != -1) MPI_Send(to_west, my_height, MPI_DOUBLE, west_rank, 1, MPI_COMM_WORLD);
        if (east_rank != -1) MPI_Recv(from_east, my_height, MPI_DOUBLE, east_rank, 1, MPI_COMM_WORLD, &status[1]);
        //south<->north boarder exchange
        if (north_rank != -1) MPI_Send(to_north, my_width, MPI_DOUBLE, north_rank, 2, MPI_COMM_WORLD);
        if (south_rank != -1) MPI_Recv(from_south, my_width, MPI_DOUBLE, south_rank, 2, MPI_COMM_WORLD, &status[2]);
        if (south_rank != -1) MPI_Send(to_south, my_width, MPI_DOUBLE, south_rank, 2, MPI_COMM_WORLD);
        if (north_rank != -1) MPI_Recv(from_north, my_width, MPI_DOUBLE, north_rank, 2, MPI_COMM_WORLD, &status[3]);
    }

    for (uint32_t i = 0; i < my_height + 2; ++i) delete[] my_block[i];
    delete[] my_block;

    if (0 == comm_rank) {
        delete[] values;
        for (int i = 0; i < height + 2; ++i) {
            delete[] values_reo[i];
        }
        delete[] values_reo;
    }

    MPI_Finalize();
    return 0;
}
