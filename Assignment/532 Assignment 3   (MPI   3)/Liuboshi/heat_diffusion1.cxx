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

void heat_update_func(int height, int width, double ** values, double ** values_next)
{
    //the border values are sources/sinks of heat, exclude halo cells
    for (uint32_t i = 1; i < height - 1; i++) {
        for (uint32_t j = 1; j < width - 1; j++) {
            values_next[i][j] = (values[i - 1][j] + values[i + 1][j] + values[i][j - 1] + values[i][j + 1]) / 4.0;
        }
    }
}

int main(int argc, char **argv) {

    if (argc != 7) {
        cout << "ERROR, incorrect arguments." << endl;
        cout << "usage:" << endl;
        cout << "\t" << argv[0] << " <simulation name : string> <height : int> <width : int> <vertical slices : int> <horizontal slices : int> <time steps : int>" << endl;
        exit(1);
    }

    string simulation_name(argv[1]);
    //the height and width of the simulation
    int height = atoi(argv[2]);
    int width = atoi(argv[3]);
    //horizontal and vertical slices will be used to determine how to divide up your MPI processes
    int vertical_slices = atoi(argv[4]);
    int horizontal_slices = atoi(argv[5]);
    //how long to run the simulation for
    int time_steps = atoi(argv[6]);

    //You will need to have a number of MPI processes equal to vertical_slices * horizontal_slices
    double **values = NULL;
    double **values_next = NULL;

    //mpi init
    MPI_Init(&argc, &argv);

    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    MPI_Status status[4];

    //slice locate
    const int row_num = comm_rank / vertical_slices;//row of this rectangle block
    const int col_num = comm_rank % vertical_slices;//col of this rectangle block
    const int row_max = vertical_slices;
    const int col_max = horizontal_slices;

    //compute neighbor ranks, invalid neighbor rank will be rank of itself
    int west_rank = row_num * row_max + (col_num - 1);
    if (0 == col_num) west_rank = -1;
    int south_rank = (row_num - 1) * row_max + col_num;
    if (0 == row_num) south_rank = -1;
    int east_rank = row_num * row_max + (col_num + 1);
    if (col_max - 1 == col_num) east_rank = -1;
    int north_rank = (row_num + 1) * row_max + col_num;
    if (row_max - 1 == row_num) north_rank = -1;

    //allocate local data and init
    const uint32_t my_height = height/vertical_slices;
    const uint32_t my_width = width/horizontal_slices;
    double ** my_block = new double*[my_height + 2];
    for (uint32_t i = 0; i < my_height + 2; ++i) {
        my_block[i] = new double[my_width + 2];
        for (uint32_t j = 0; j < my_width + 2; ++j) {
            my_block[i][j] = 0.0;
        }
    }

    //put heat source
    for (int i = 0; i < my_height; ++i) {
        from_west[i] = (0 == col_num) ? 1.0 : 0.0;
    }
    for (int i = 0; i < my_height; ++i) {
        from_west[i] = (col_max - 1 == col_num) ? -1.0 : 0.0;
    }
    for (int i = 0; i < my_width; ++i) {

    }
    if (0 == col_num) {
        for (int i = 0; i < my_height; ++i) {
            from_west[i] = 1.0;
        }
    } else if (col_max - 1 == col_num) {
        for (int i = 0; i < my_height; ++i) {
            from_east[i] = -1.0;
        }
    }
    if (0 = row_num) {
        for (int i = 0; i < my_width; ++i) {
            from_south[i] = 0.0;
        }
    } else if (col_max - 1 == col_num) {
        for (int i = 0; i < my_width; ++i) {
            from_north[i] = 0.0;
        }
    }

    //allocate result array
    double my_result[my_height*my_width];
    for (uint32_t i = 0; i < my_height*my_width; ++i) {
        my_result[i] = 0.0;
    }

    //allocate halo array for transfer
    double from_east[my_height];
    double to_east[my_height];
    double from_west[my_height];
    double to_west[my_height];
    double from_north[my_width];
    double to_north[my_width];
    double from_south[my_width];
    double to_south[my_width];

    for (int steps = 0; steps < time_steps; steps++)
    {
        //fill in halo
        for (int i = 0; i < my_height; ++i) {
            from_west[i] = 1.0;
        }
        for (int i = 0; i < my_height; ++i) {
            from_east[i] = -1.0;
        }
        for (int i = 0; i < my_width; ++i) {
            from_south[i] = 0.0;
        }
        for (int i = 0; i < my_width; ++i) {
            from_north[i] = 0.0;
        }


        //perform a step to simulate the heat diffusion
        for (uint32_t i = 0; i < my_height; ++i) {
            for (uint32_t j = 0; j < my_width; ++j) {
                my_result[i*my_width + j] = (my_block[i][j + 1] + values[i + 2][j + 1] + values[i + 1][j] + values[i + 1][j + 2]) / 4.0;
            }
        }

        //fill in halo to transfer
        for (int i = 0; i < my_height; ++i) {
            to_east[i] = my_result[i * my_width + (my_width - 1)];
        }
        for (int i = 0; i < my_height; ++i) {
            to_west[i] = my_result[i*my_width + 0];
        }
        for (int i = 0; i < my_width; ++i) {
            to_north[i] = my_result[(my_height - 1) * my_width + i];
        }
        for (int i = 0; i < my_width; ++i) {
            to_south[i] = my_result[0 + i];
        }

        //halo tansfer via MPI
        //east<->west boarder exchange
        MPI_Send(to_east, my_height, MPI_DOUBLE, east_rank, 1, MPI_COMM_WORLD);
        MPI_Recv(from_west, my_height, MPI_DOUBLE, west_rank, 1, MPI_COMM_WORLD, &status[0]);
        MPI_Send(to_west, my_height, MPI_DOUBLE, west_rank, 1, MPI_COMM_WORLD);
        MPI_Recv(from_east, my_height, MPI_DOUBLE, east_rank, 1, MPI_COMM_WORLD, &status[1]);
        //south<->north boarder exchange
        MPI_Send(to_north, my_width, MPI_DOUBLE, north_rank, 2, MPI_COMM_WORLD);
        MPI_Recv(from_south, my_width, MPI_DOUBLE, south_rank, 2, MPI_COMM_WORLD, &status[2]);
        MPI_Send(to_south, my_width, MPI_DOUBLE, south_rank, 2, MPI_COMM_WORLD);
        MPI_Recv(from_north, my_width, MPI_DOUBLE, north_rank, 2, MPI_COMM_WORLD, &status[3]);


    }
}
