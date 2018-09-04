#include <cstdlib>
#include <iostream>
#include <ostream>
#include <sstream>  //for stringstream/ostringstream/istringstream

#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "unistd.h"

using std::cout;
using std::endl;
using std::ostringstream;

int main(int argc, char **argv) {
    int     comm_sz;    /* number of processes  */
    int     my_rank;    /* process rank         */

    //Initialize the command line arguments on every process
    MPI_Init(&argc, &argv);

    //Get the number of processes in MPI_COMM_WORLD, and put 
    //it in the 'comm_sz" variable; ie., how many processes
    //are running this program (the same as what you put in the 
    //-np argument).
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    //Get the rank of this particular process in MPI_COMM_WORLD,
    //and put it in the 'my_rank' variable -- ie., what number
    //is this process
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int token[1];

    if (my_rank == 0) {
        int terminate_tags_received = 0;
        int message_count = 0;

        while (true) {
            MPI_Request request;
            MPI_Irecv(&token, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);

            cout << "[RANK " << my_rank << "] set up asynchronous recv!" << endl;

            MPI_Status status;
            MPI_Wait(&request, &status);
            message_count++;

            cout << "[RANK " << my_rank << "] received a message from process " << status.MPI_SOURCE << " with tag " << status.MPI_TAG << " and value: " << token[0] << ", total messages received: " << message_count << endl;

            if (status.MPI_TAG == 1) terminate_tags_received++;

            if (terminate_tags_received == (comm_sz - 1)) break;
        }

    } else {
        int max_messages = 5;
        for (int i = 0; i < max_messages; i++) {
            token[0] = drand48() * 100;
            int tag = 0;
            if (i == (max_messages - 1)) tag = 1;

            MPI_Request request;
            MPI_Isend(&token, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);

            cout << "[RANK " << my_rank << "] set up asynchronous send!" << endl;

            MPI_Status status;
            MPI_Wait(&request, &status);

            cout << "[RANK " << my_rank << "] sent message!" << endl;

            double sleep_time = 5 * drand48();

            cout << "[RANK " << my_rank << "] sleeping for " << sleep_time << " seconds!" << endl;
            sleep(sleep_time);
        }
    }

    MPI_Finalize();
}