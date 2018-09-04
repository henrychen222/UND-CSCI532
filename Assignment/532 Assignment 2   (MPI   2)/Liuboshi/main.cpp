#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#if (__cplusplus == 201103L)
#include <unordered_map>
#else
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#endif

using namespace std;

#include <mpi/mpi.h>
#include <string.h>

int main(int argc, char** argv)
{
    int comm_sz;
    int comm_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    MPI_Request comm_request[comm_sz];
    MPI_Status comm_status[comm_sz];

    //read seed file in master and dispatch to all processes
    unordered_map<string, int> seeds;
    char * my_seed_list = nullptr;
    int my_seed_num = 0;
    int all_seed_num = 0;

    int proc_seed_counts[comm_sz];
    memset(proc_seed_counts, 0, sizeof(int)*comm_sz);
    string record;
    vector<string> chromosome;
    vector<string> chr_tag;
    unordered_map<int, int> * result;
    int stop_command[comm_sz] = {0};

    //step 1 read files, master read seed file and slave process read chr file
    string read_line;
    if(0 == comm_rank)
    {
        //master procees read the seed files, store seed into unordered_map(merge the homo-sequences) and serialize to an array
        seeds.rehash(10000000);
        int amount = 0;
        //read seed file in the master
        for (int idx = 0; idx < argc; idx++)
        {
            if (0 == strcmp(argv[idx], "--seed_files"))
            {
                //read each file under the --seed_files
                while (++idx < argc)
                {
                    //detected another --, finish reading files
                    if (strlen(argv[idx]) > 2 && argv[idx][0] == '-' && argv[idx][1] == '-') break;
                    ifstream in(argv[idx]);
                    if (in.is_open())
                    {
                        //read every line and find out the seed sequence
                        string last_read("");
                        while (getline(in, read_line))
                        {

                            if(1 == read_line.size() && '+' == read_line.at(0))
                            {
                                //save the seed that have 40 chars without a 'N' before a line "+"
                                if (40 == last_read.size() && string::npos == last_read.find('N'))
                                {
                                    seeds[last_read] += 1;
                                    amount++;
                                }
                            }
                            last_read = read_line;
                        }
                    }
                    in.close();
                }
                //quit the loop if all seed file have been read
                break;
            }
        }
        //print logs
        all_seed_num = seeds.size();
        cout << "In process[" << comm_rank << "] seeds files have " << all_seed_num << " seeds " << amount << endl;
    }
    //all processes read chromesome files
    for (int idx = 0; idx < argc; idx++)
    {
        if (0 == strcmp(argv[idx], "--chr_files"))
        {
            while (++idx < argc)
            {
                if (strlen(argv[idx]) > 2 && argv[idx][0] == '-' && argv[idx][1] == '-') break;
                ifstream in(argv[idx]);
                if (in.is_open())
                {
                    //record chr id
                    getline(in, read_line);
                    chr_tag.push_back(read_line.substr(1));
                    //read in and convert
                    while (getline(in, read_line))
                    {
                        //all to upper case
                        for (auto & var : read_line)
                        {
                            var = toupper(var);
                        }
                        //append
                        record.append(read_line);
                    }
                }
                in.close();
                chromosome.push_back(record);
                record.clear();
                //print logs
                cout << "In process[" << comm_rank << "] has read " << chr_tag.size() << " chr files! " << endl;
            }
        }
    }
    result = new unordered_map<int,int>[chr_tag.size()];
    for(int idx = 0; idx < chr_tag.size(); ++idx)
    {
        result[idx].clear();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&all_seed_num, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double count = 0.0;
    int prev_count = 0;
    proc_seed_counts[0] = 0;
    for (int idx = 1; idx < comm_sz; ++idx)
    {
        prev_count = (int)count;
        count += (double)all_seed_num/(double)(comm_sz - 1);
        proc_seed_counts[idx] = (int)count - prev_count;
    }

    //alloc memory space for each process and dispatch seeds to each process
    if(0 == comm_rank)
    {
        //master extract all seeds then send to the slave processes
        my_seed_num = all_seed_num;
        my_seed_list = new char[all_seed_num * 40];
        int idx = 0;
        for(const auto & var : seeds)
        {
            memcpy(my_seed_list + idx * 40, var.first.c_str(), sizeof(char) * 40);
            idx++;
        }
        int start_pos = 0;
        for (int idx = 1; idx < comm_sz; ++idx)
        {
            MPI_Isend(my_seed_list + start_pos, proc_seed_counts[idx] * 40, MPI_CHAR, idx, 100, MPI_COMM_WORLD, &comm_request[0]);
            MPI_Wait(&comm_request[0], &comm_status[0]);
            start_pos += proc_seed_counts[idx] * 40;
        }
    }
    else
    {
        //slave processes receive seeds
        my_seed_num = proc_seed_counts[comm_rank];
        my_seed_list = new char[my_seed_num * 40];
        memset(my_seed_list, '\0', sizeof(char)*my_seed_num * 40);
        MPI_Irecv(my_seed_list, my_seed_num * 40, MPI_CHAR, 0, 100, MPI_COMM_WORLD, &comm_request[comm_rank]);
        MPI_Wait(&comm_request[comm_rank], &comm_status[comm_rank]);

        string key_stream(my_seed_list);
        //put the seed into the unordered_map
        for (int idx = 0; idx < my_seed_num; ++idx)
        {
            seeds[key_stream.substr(idx * 40, 40)] = 1;
        }

        std::cout << "Proc " << comm_rank << " get " << my_seed_num << " seeds!  ";
        cout << my_seed_list[0] << my_seed_list[1] << " ... " << my_seed_list[my_seed_num*40 - 2] << my_seed_list[my_seed_num*40 - 1] << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    const int comm_buf_len = 40 + sizeof(size_t);
    char t_r_buff[comm_buf_len];

    if (0 == comm_rank)
    {
        //master procee is in charge of handling the message from each slave process and write file
        //loop to receive message
        while(true)
        {
            //receive a message from any source and any tag
            MPI_Irecv(t_r_buff, comm_buf_len, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &comm_request[0]);
            MPI_Wait(&comm_request[0], &comm_status[0]);
            //if this message is a match message then tag should greater then 0 or 0
            if (comm_status[0].MPI_TAG < 100)
            {
                size_t * val = reinterpret_cast<size_t *>(t_r_buff + 40 * sizeof(char));
                string lookup(t_r_buff);
                lookup = lookup.substr(0, 40);
                result[comm_status[0].MPI_TAG][(int)(*val)] += seeds[lookup];
            }
            if (101 == comm_status[0].MPI_TAG)
            {
                stop_command[comm_status[0].MPI_SOURCE] = 1;
            }
            int sum = 0;
            for (int idx = 0; idx < comm_sz; ++idx)
            {
                sum += stop_command[idx];
            }
            if (sum >= comm_sz - 1)
            {
                break;
            }
        }
        ofstream out("./output.txt", ios_base::out | ios_base::trunc);
        for (int i = 0; i < chr_tag.size(); ++i)
        {
            for(const auto & var : result[i])
            {
                out << chr_tag.at(i) << ',' << var.first << ',' << var.second << endl;
            }
        }
        out.close();
    }
    else
    {
        //slave processes start to process and master process handle the receive
        //for each chromosome, run a compare procedure
        size_t send_count = 0;
        for(int idx = 0; idx < chr_tag.size(); idx++)
        {
            string & seq = chromosome.at(idx);
            string test_str;
            size_t cur = 0;
            while (cur < seq.size() - 40)
            {
                //test 'N', if there exist N, adjust the cur to the next pos of last N.
                test_str = seq.substr(cur, 40);
                size_t pos = test_str.find_last_of('N');
                if(pos != string::npos)
                {
                    cur += pos;
                    cur++;
                }
                else
                {
                    //No N exist, search in unordered_map, report when matched
                    if(seeds.end() != seeds.find(test_str))
                    {
                        memcpy(t_r_buff, test_str.c_str(), sizeof(char)*40);
                        memcpy(t_r_buff + 40 * sizeof(char), &cur, sizeof(size_t));
                        //tag indicate which chr matched this seed, the tag means the chr id which is >=0
                        //cout << "\tProc[" << comm_rank << "] matched at " << cur << "/" << seq.size() << " " << test_str << endl;
                        MPI_Isend(t_r_buff, comm_buf_len, MPI_CHAR, 0, idx, MPI_COMM_WORLD, &comm_request[comm_rank]);
                        MPI_Wait(&comm_request[comm_rank], &comm_status[comm_rank]);
                        send_count++;
                    }
                    cur++;
                }
            }
        }
        //after all chr files complete matching, send stop
        MPI_Isend(t_r_buff, comm_buf_len, MPI_CHAR, 0, 101, MPI_COMM_WORLD, &comm_request[comm_rank]);
        cout << "Proc[" << comm_rank << "] sent stop code 100 to MASTER after (" << send_count << ") matches sent!" << endl;
    }

    delete[] my_seed_list;
    delete[] result;

    //finalize
    MPI_Finalize();

    return 0;
}
