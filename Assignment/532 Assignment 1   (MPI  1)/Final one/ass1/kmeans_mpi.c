#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>

#include <mpi/mpi.h>

double dist_2(double x, double y, double z, double a, double b, double c)
{
    return (x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c);
}

int main(int argc, char** argv)
{
    int comm_sz;
    int comm_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    //get the cluster num
    int centroid_num = 0;            //number of cluster
    double * centroids = NULL;       //centroid data
    double * centroid_avg = NULL;
    for(int idx = 1; idx < argc; ++idx)
    {
        if(0 == strcmp(argv[idx], "--cluster_num"))
        {
            centroid_num = atoi(argv[++idx]);
            centroids = (double *)malloc(sizeof(double)*centroid_num*3);
            memset(centroids, 0, sizeof(double)*centroid_num*3);
        }
    }

    int all_star_num = 0;                //number of stars
    double * all_stars_xyz = NULL;       //the xyz star coordination data
    double * all_stars_lbr = NULL;       //the lbr coordination orignal data
    int * all_assignment = NULL;         //identify every star with a cluster number

    //master process in charge of read file and generate random centroids
    if(0 == comm_rank)
    {
        //read star files and get star nums
        for (int idx = 0; idx < argc; idx++)
        {
            if(0 == strcmp(argv[idx], "--star_files"))
            {
                while(++idx < argc)
                {
                    if(strlen(argv[idx]) > 2 && argv[idx][0] == '-' && argv[idx][1] == '-') break;
                    FILE * fp = fopen(argv[idx], "r");
                    int digits = 0;
                    if(NULL != fp)
                    {
                        fscanf(fp, "%d", &digits);
                        all_star_num += digits;
                    }
                    fclose(fp);
                }
                break;
            }
        }
        //allocate memory for star data and reslut data
        all_stars_xyz = (double *)malloc(sizeof(double)*all_star_num*3);
        all_stars_lbr = (double *)malloc(sizeof(double)*all_star_num*3);
        all_assignment = (int *)malloc(sizeof(int)*all_star_num);
        centroid_avg = (double *)malloc(sizeof(double)*centroid_num*4);

        memset(all_stars_xyz, 0, sizeof(double)*all_star_num*3);
        memset(all_stars_lbr, 0, sizeof(double)*all_star_num*3);
        memset(all_assignment, 0, sizeof(int)*all_star_num);
        memset(centroid_avg, 0, sizeof(double)*centroid_num*4);
    }

    //broadcast star number
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&all_star_num, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(0 == comm_rank)
    {
        //read and analysis star files
        double max_x = -DBL_MIN, min_x = DBL_MAX;
        double max_y = -DBL_MIN, min_y = DBL_MAX;
        double max_z = -DBL_MIN, min_z = DBL_MAX;
        int occu = 0;
        int read_num = 0;
        double l = 0.0, b = 0.0, r = 0.0;
        double xx = 0.0, yy = 0.0, zz = 0.0;
        //read star files and load star data
        for(int idx = 1; idx < argc; ++idx)
        {
            if(0 == strcmp(argv[idx], "--star_files"))
            {
                ++idx;
                while(idx < argc)
                {
                    if(strlen(argv[idx]) > 2 && argv[idx][0] == '-' && argv[idx][1] == '-') break;
                    //open star file
                    FILE * fp = fopen(argv[idx], "r");
                    //ommit bad files
                    if(NULL == fp) continue;
                    //fetch and transform star datas
                    fscanf(fp, "%d", &occu);
                    while(EOF != fscanf(fp, "%lf %lf %lf", &l, &b, &r))
                    {
                        //record input stars
                        all_stars_lbr[read_num * 3 + 0] = l;
                        all_stars_lbr[read_num * 3 + 1] = b;
                        all_stars_lbr[read_num * 3 + 2] = r;
                        //convert lbr (galactic) to x y z (cartesian) [formula from star_visualizier]
                        l = l * 3.14159265358979323846 / 180;
                        b = b * 3.14159265358979323846 / 180;
                        xx = r * cos(b) * sin(l);
                        yy = r * cos(l) * cos(b);
                        zz = r * sin(b);
                        all_stars_xyz[read_num * 3 + 0] = xx;
                        all_stars_xyz[read_num * 3 + 1] = yy;
                        all_stars_xyz[read_num * 3 + 2] = zz;
                        //do statistics
                        read_num++;
                        min_x = xx < min_x ? xx : min_x; max_x = xx > max_x ? xx : max_x;
                        min_y = yy < min_y ? yy : min_y; max_y = yy > max_y ? yy : max_y;
                        min_z = zz < min_z ? zz : min_z; max_z = zz > max_z ? zz : max_z;
                    }
                    fclose(fp);
                    ++idx;
                }
            }
        }

        //print logs
        printf("Master>>Stars succesfully loaded.\n");
        printf("Master>>\tstars exist: %d\n", all_star_num);
        printf("Master>>\tstars read: %d\n", read_num);
        printf("Master>>\tmin / max x val: %lf / %lf\n", min_x, max_x);
        printf("Master>>\tmin / max y val: %lf / %lf\n", min_y, max_y);
        printf("Master>>\tmin / max z val: %lf / %lf\n", min_z, max_z);

        //generate random (x y z) for kmeans algorithm
        for(int cnt = 0; cnt < centroid_num; ++cnt)
        {
            srand(cnt*100+13);
            centroids[cnt*3+0] = (max_x - min_x)*(rand()%10/10.0) + min_x;
            srand(cnt*200+31);
            centroids[cnt*3+1] = (max_y - min_y)*(rand()%10/10.0) + min_y;
            srand(cnt*250+47);
            centroids[cnt*3+2] = (max_z - min_z)*(rand()%10/10.0) + min_z;
        }

        printf("K-Means Centroid Init:\n");
        for(int cnt = 0; cnt < centroid_num; ++cnt)
        {
            printf("c%d(%.5lf,%.5lf,%.5lf) ", cnt, centroids[cnt*3+0], centroids[cnt*3+1], centroids[cnt*3+2]);
        }
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(centroids, centroid_num*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast(all_stars_xyz, all_star_num*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printf("P[%d/%d] Bcast centroid ", comm_rank, comm_sz);
    for(int cnt = 0; cnt < centroid_num; ++cnt)
    {
        printf("c%d(%.5lf,%.5lf,%.5lf) ", cnt, centroids[cnt*3+0], centroids[cnt*3+1], centroids[cnt*3+2]);
    }
    printf("\n");
    MPI_Barrier(MPI_COMM_WORLD);

    //now start to slice stars and assignments into different slices
    int * star_counts = (int *)malloc(sizeof(int)*comm_sz);
    int * star_displs = (int *)malloc(sizeof(int)*comm_sz);
    int * assi_counts = (int *)malloc(sizeof(int)*comm_sz);
    int * assi_displs = (int *)malloc(sizeof(int)*comm_sz);
    double count = 0.0;
    int prev_count = 0;
    for(int idx = 0; idx < comm_sz; ++idx)
    {
        //slice assignments
        assi_displs[idx] = (int)count;
        prev_count = (int)count;
        count += (double)all_star_num/(double)comm_sz;
        assi_counts[idx] = (int)count - (int)prev_count;
        //slice stars
        star_counts[idx] = 3 * assi_counts[idx];
        star_displs[idx] = 3 * assi_displs[idx];
    }

    //alloc memory space for each process
    double * my_stars = (double *)malloc(sizeof(double)*star_counts[comm_rank]);
    memset(my_stars, 0, sizeof(double)*star_counts[comm_rank]);
    int * my_assis = (int *)malloc(sizeof(int)*assi_counts[comm_rank]);
    memset(my_assis, 0, sizeof(int)*assi_counts[comm_rank]);
    double sse = 0.0;

    //scatter stars
    MPI_Barrier(MPI_COMM_WORLD);

    printf("P[%d/%d] assi(star)slice [", comm_rank, comm_sz);
    for(int idx = 0; idx < comm_sz; ++idx)
    {
        printf("%d(%d)  ", star_counts[idx], assi_counts[idx]);
    }
    printf("]\n");

    printf("P[%d/%d] displacese [", comm_rank, comm_sz);
    for(int idx = 0; idx < comm_sz; ++idx)
    {
        printf("%d(%d)  ", star_displs[idx], assi_displs[idx]);
    }
    printf("]\n");

    MPI_Scatterv(all_stars_xyz, star_counts, star_displs, MPI_DOUBLE, my_stars, star_counts[comm_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("P[%d/%d] Scatter %d stars\n", comm_rank, comm_sz, star_counts[comm_rank]);

    MPI_Barrier(MPI_COMM_WORLD);

    int loop_count = 0;
    int loop_flag = 1;

    while(1 == loop_flag)
    {
        //initial clustering: inite assi and sse (select the nearest centroid as each star's cluster number).
        for(int s_id = 0; s_id < star_counts[comm_rank] / 3; ++s_id)
        {
            double min_dist_2 = DBL_MAX;
            for(int c_id = 0; c_id < centroid_num; ++c_id)
            {
                //calulate dist
                double my_dist_2 = dist_2(centroids[c_id*3+0], centroids[c_id*3+1], centroids[c_id*3+2], my_stars[s_id*3+0], my_stars[s_id*3+1], my_stars[s_id*3+2]);
                //if star tested is closer to this centroid, update its cluster number stored in assi
                if(my_dist_2 < min_dist_2)
                {
                    //update minima dist
                    min_dist_2 = my_dist_2;
                    //update the assignment array
                    my_assis[s_id] = c_id;
                }
            }
        }

        //gather assignment
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gatherv(my_assis, assi_counts[comm_rank], MPI_INT, all_assignment, assi_counts, assi_displs, MPI_INT, 0, MPI_COMM_WORLD);

        if(0 == comm_rank)
        {
            //inc the loop count
            loop_count++;

            //calcu the new sse
            double new_sse = 0.0;
            memset(centroid_avg, 0, sizeof(double)*centroid_num*4);
            for(int cnt = 0; cnt < all_star_num; ++cnt)
            {
                new_sse += dist_2(all_stars_xyz[cnt*3+0], all_stars_xyz[cnt*3+1], all_stars_xyz[cnt*3+2], centroids[all_assignment[cnt]*3+0], centroids[all_assignment[cnt]*3+1], centroids[all_assignment[cnt]*3+2]);
                centroid_avg[all_assignment[cnt]*4+0] += all_stars_xyz[cnt*3+0];
                centroid_avg[all_assignment[cnt]*4+1] += all_stars_xyz[cnt*3+1];
                centroid_avg[all_assignment[cnt]*4+2] += all_stars_xyz[cnt*3+2];
                centroid_avg[all_assignment[cnt]*4+3] += 1.0;
            }
            double conv = fabs(new_sse - sse);
            sse = new_sse;

            //update centroid based on initial cluster result;
            for(int cnt = 0; cnt < centroid_num; ++cnt)
            {
                if(0 != centroid_avg[cnt*4+3])
                {
                    centroids[cnt*3+0] = centroid_avg[cnt*4+0] / centroid_avg[cnt*4+3];
                    centroids[cnt*3+1] = centroid_avg[cnt*4+1] / centroid_avg[cnt*4+3];
                    centroids[cnt*3+2] = centroid_avg[cnt*4+2] / centroid_avg[cnt*4+3];
                }
            }

            printf("\nK-Means Centroid Now(%dth iter) delta_sse = %lf:\n", loop_count, conv);
            for(int cnt = 0; cnt < centroid_num; ++cnt)
            {
                printf("  c%d(%.5lf,%.5lf,%.5lf) ", cnt, centroids[cnt*3+0], centroids[cnt*3+1], centroids[cnt*3+2]);
            }
            printf("\n");

            if(conv < 0.001)
            {
                loop_flag = 0;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(centroids, centroid_num*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&loop_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        loop_count++;
    }

    if(0 == comm_rank)
    {
        //printing message and logging
        printf("\nCluster assignmet:\n");
        for(int cnt = 0; cnt < centroid_num; ++cnt)
        {
            printf("\tc%d[%.5lf,%.5lf,%.5lf] with %d in %d stars %.2lf%%\n", cnt, centroids[cnt*3+0], centroids[cnt*3+1], centroids[cnt*3+2], (int)centroid_avg[cnt*4+3], all_star_num, centroid_avg[cnt*4+3]*100.0/all_star_num);
        }

        //write to files
        //write non-mpi kmeans result to files
        system("rm -rf ./clusters_mpi");
        mkdir("./clusters_mpi", 0777);
        FILE ** fa = (FILE **)malloc(sizeof(FILE*)*centroid_num);
        for(int cnt = 0; cnt < centroid_num; ++cnt)
        {
            fa[cnt] = NULL;
        }
        char fn[128] = {'\0'};
        //create and open files
        for(int cnt = 0; cnt < centroid_num; ++cnt)
        {
            if(0 != (int)centroid_avg[cnt*4+3])
            {
                sprintf(fn, "./clusters_mpi/stars-%d.txt", cnt);
                fa[cnt] = fopen(fn, "w");
                fprintf(fa[cnt], "%d\n", (int)centroid_avg[cnt*4+3]);
            }
        }
        //write data of each cluster, Ommit the empty file
        for(int cnt = 0; cnt < all_star_num; ++cnt)
        {
            if(NULL != fa[all_assignment[cnt]])
            {
                fprintf(fa[all_assignment[cnt]], "%f %f %f\n", all_stars_lbr[cnt*3+0], all_stars_lbr[cnt*3+1], all_stars_lbr[cnt*3+2]);
            }
        }
        //close files
        for(int cnt = 0; cnt < centroid_num; cnt++)
        {
            if(NULL != fa[cnt])
            {
                fclose(fa[cnt]);
            }
        }
        free(fa);

    }

    //dellocate memory
    if(NULL != all_stars_xyz) { free(all_stars_xyz); all_stars_xyz = NULL; }
    if(NULL != all_stars_lbr) { free(all_stars_lbr); all_stars_lbr = NULL; }
    if(NULL != centroids) { free(centroids); centroids = NULL; }
    if(NULL != centroid_avg) { free(centroid_avg); centroid_avg = NULL; }
    if(NULL != all_assignment) { free(all_assignment); all_assignment = NULL; }
    if(NULL != star_counts) { free(star_counts); star_counts = NULL; }
    if(NULL != star_displs) { free(star_displs); star_displs = NULL; }
    if(NULL != assi_counts) { free(assi_counts); assi_counts = NULL; }
    if(NULL != assi_displs) { free(assi_displs); assi_displs = NULL; }
    if(NULL != my_stars) { free(my_stars); my_stars = NULL; }
    if(NULL != my_assis) { free(my_assis); my_assis = NULL; }

    //finalize
    MPI_Finalize();

    return 0;
}




