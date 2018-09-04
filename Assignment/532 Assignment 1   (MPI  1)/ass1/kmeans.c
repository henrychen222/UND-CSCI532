#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>

//for atgument reading
static int centroid_num = 0;            //number of clusters
static int star_num = 0;                //number of stars
static int star_read_num = 0;           //number read from star file
static int star_file_num = 0;           //number of star file

//for initialization
static double max_x = -DBL_MIN;
static double min_x = DBL_MAX;
static double max_y = -DBL_MIN;
static double min_y = DBL_MAX;
static double max_z = -DBL_MIN;
static double min_z = DBL_MAX;

//for cluster data record
static double * stars = NULL;           //the xyz star coordination data
static double * records = NULL;         //the lbr coordination orignal data
static int * cluster = NULL;            //identify every star with a cluster number
static int * cluster_count = NULL;      //count the number in each cluster

//for kmeans iteration
static double * centroid_x = NULL;
static double * centroid_y = NULL;
static double * centroid_z = NULL;
static double * centroid_x_itr = NULL;
static double * centroid_y_itr = NULL;
static double * centroid_z_itr = NULL;
static int * centroid_n_itr = NULL;

double dist_2(double x, double y, double z, double a, double b, double c)
{
    return (x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c);
}

//Lloyd's algorithm
int lloyd_kmeans(double stars[], int star_num, int c_num, double c_x[], double c_y[], double c_z[], int result[])
{
    memset(result, 0, sizeof(int)*star_num);
    int loop_counter = 1;
    double sse = 0.0;
    for(int idx = 0; idx < c_num; ++idx)
    {
        centroid_x_itr[idx] = 0.0;
        centroid_y_itr[idx] = 0.0;
        centroid_z_itr[idx] = 0.0;
        centroid_n_itr[idx] = 0;
    }
    struct timeval tik, toc;

    //tik
    gettimeofday(&tik, NULL);

    //initial try to cluster: find initial sse, cluster id(select the nearest centroid as star's cluster_id).
    for(int star_id = 0; star_id < star_num; ++star_id)
    {
        double min_dist_2 = DBL_MAX;
        for(int c_id = 0; c_id < c_num; ++c_id)
        {
            //calulate dist
            double my_dist_2 = dist_2(c_x[c_id], c_y[c_id], c_z[c_id], stars[star_id*3+0], stars[star_id*3+1], stars[star_id*3+2]);
            //if current star tested to be more close to this centroid
            if(my_dist_2 < min_dist_2)
            {
                //update minima dist
                min_dist_2 = my_dist_2;
                //update the cluster
                result[star_id] = c_id;
            }
        }
        centroid_x_itr[result[star_id]] += stars[star_id * 3 + 0];
        centroid_y_itr[result[star_id]] += stars[star_id * 3 + 1];
        centroid_z_itr[result[star_id]] += stars[star_id * 3 + 2];
        centroid_n_itr[result[star_id]] += 1;
    }

    //update centroid based on initial cluster result;
    for(int cnt = 0; cnt < c_num; ++cnt)
    {
        if(0 != centroid_n_itr[cnt])
        {
            c_x[cnt] = centroid_x_itr[cnt] / centroid_n_itr[cnt];
            c_y[cnt] = centroid_y_itr[cnt] / centroid_n_itr[cnt];
            c_z[cnt] = centroid_z_itr[cnt] / centroid_n_itr[cnt];
        }
    }
    //calculate sse
    for(int idx = 0; idx < star_num; ++idx)
    {
        sse += dist_2(stars[idx*3+0], stars[idx*3+1], stars[idx*3+2], c_x[result[idx]], c_y[result[idx]], c_z[result[idx]]);
    }

    //loop and update centroid and sse until convergence
    while(1)
    {
        //update cluster id for every star.
        //init statistic vars
        for(int idx = 0; idx < c_num; ++idx)
        {
            centroid_x_itr[idx] = 0.0;
            centroid_y_itr[idx] = 0.0;
            centroid_z_itr[idx] = 0.0;
            centroid_n_itr[idx] = 0;
        }

        //select the nearest centroid as star's cluster_id.(on init configuration, all star belongs to cluster 1)
        for(int star_id = 0; star_id < star_num; ++star_id)
        {
            double min_dist_2 = DBL_MAX;
            for(int c_id = 0; c_id < c_num; ++c_id)
            {
                double my_dist_2 = dist_2(c_x[c_id], c_y[c_id], c_z[c_id], stars[star_id*3+0], stars[star_id*3+1], stars[star_id*3+2]);
                //if current star tested to be more close to
                if(my_dist_2 < min_dist_2)
                {
                    //update minima dist
                    min_dist_2 = my_dist_2;
                    //update the cluster
                    result[star_id] = c_id;
                }
            }
            centroid_x_itr[result[star_id]] += stars[star_id * 3 + 0];
            centroid_y_itr[result[star_id]] += stars[star_id * 3 + 1];
            centroid_z_itr[result[star_id]] += stars[star_id * 3 + 2];
            centroid_n_itr[result[star_id]] += 1;
        }

        //update centroid based on new cluster result;
        for(int cnt = 0; cnt < c_num; ++cnt)
        {
            if(0 != centroid_n_itr[cnt])
            {
                c_x[cnt] = centroid_x_itr[cnt] / centroid_n_itr[cnt];
                c_y[cnt] = centroid_y_itr[cnt] / centroid_n_itr[cnt];
                c_z[cnt] = centroid_z_itr[cnt] / centroid_n_itr[cnt];
            }
        }

        //calculate new sse and update count of each cluster
        double new_sse = 0.0;
        memset(cluster_count, 0, sizeof(int)*c_num);
        for(int idx = 0; idx < star_num; ++idx)
        {
            new_sse += dist_2(stars[idx*3+0], stars[idx*3+1], stars[idx*3+2], c_x[result[idx]], c_y[result[idx]], c_z[result[idx]]);
            ++cluster_count[result[idx]];
        }
        double delta_sse = fabs(new_sse - sse);

        printf("\nK-Means Centroid Now(%dth iter) delta_sse = %lf:\n", loop_counter, delta_sse);
        for(int cnt = 0; cnt < c_num; ++cnt)
        {
            printf("\tc%d(%.5lf,%.5lf,%.5lf)  ", cnt, c_x[cnt], c_y[cnt], c_z[cnt]);
        }

        loop_counter++;

        //break if convergence or update sse
        if(delta_sse <= 0.001)
        {
            break;
        }
        else
        {
            sse = new_sse;
        }
    }

    //toc when convergence
    gettimeofday(&toc, NULL);

    //printing message and logging
    printf("\nCluster result:\n");
    for(int cnt = 0; cnt < c_num; ++cnt)
    {
        printf("\tc%d[%.5lf,%.5lf,%.5lf] with %d in %d stars %.2lf%%\n", cnt, c_x[cnt], c_y[cnt], c_z[cnt], cluster_count[cnt], star_num, cluster_count[cnt]*100.0/star_num);
    }
    printf("\nTime elapsed: %10.6lfs\n", (toc.tv_sec - tik.tv_sec) + (toc.tv_usec - tik.tv_usec)/1000000.0);

    return 0;
}

int main(int argc, char** argv)
{
    //step-0 parse argument and do some init
    if(argc > 2)
    {
        //handle cluster_num
        centroid_num = atoi(argv[1]);

        //handle star files
        for(int idx = 2; idx < argc; ++idx)
        {
            if(strlen(argv[idx]) > 0)
            {
                FILE * fp = fopen(argv[idx], "r");
                int digits = 0;
                if(NULL != fp)
                {
                    fscanf(fp, "%d", &digits);
                    star_num += digits;
                    star_file_num++;
                }
                fclose(fp);
            }
        }

        //print errors or messages
        printf("Arguments succesfully parsed.\n");
        printf("\tstar file number: %d\n", star_file_num);
        printf("\tstar number: %d\n", star_num);
        printf("\tcluster number: %d\n", centroid_num);
        if(centroid_num < 1)
        {
            printf("\tERROR: Strange cluster num [%d]... ABORT\n", centroid_num);
            return 1;
        }
        if(0 == star_file_num)
        {
            printf("\tERROR: No valid star files... ABORT\n");
            return 1;
        }
        if(0 == star_num)
        {
            printf("\tERROR: No stars to be clustered... ABORT\n");
            return 1;
        }
    }
    else
    {
        printf("Wrong arguments... ABORT!\n");
        return 1;
    }

    //step-1 Allocate memory space for star data
    stars = (double *)malloc(sizeof(double)*star_num*3);
    records = (double *)malloc(sizeof(double)*star_num*3);
    centroid_x = (double *)malloc(sizeof(double)*centroid_num);
    centroid_y = (double *)malloc(sizeof(double)*centroid_num);
    centroid_z = (double *)malloc(sizeof(double)*centroid_num);
    centroid_x_itr = (double *)malloc(sizeof(double)*centroid_num);
    centroid_y_itr = (double *)malloc(sizeof(double)*centroid_num);
    centroid_z_itr = (double *)malloc(sizeof(double)*centroid_num);
    centroid_n_itr = (int *)malloc(sizeof(int)*centroid_num);
    cluster = (int *)malloc(sizeof(int)*star_num);
    cluster_count = (int *)malloc(sizeof(int)*centroid_num);

    //step-2 read star files
    int occu = 0;
    double l = 0.0, b = 0.0, r = 0.0;
    double xx = 0.0, yy = 0.0, zz = 0.0;
    //open every star file to read
    for(int idx = 2; idx < argc; ++idx)
    {
        if(strlen(argv[idx]) > 0)
        {
            //open star file
            FILE * fp = fopen(argv[idx], "r");
            //ommit bad files
            if(NULL == fp) continue;
            //fetch and transform star datas
            fscanf(fp, "%d", &occu);
            while(EOF != fscanf(fp, "%lf %lf %lf", &l, &b, &r))
            {
                //convert degrees to radians
                l = l * 3.14159265358979323846 / 180;
                b = b * 3.14159265358979323846 / 180;
                //convert lbr (galactic) to x y z (cartesian)
                xx = r * cos(b) * sin(l);
                yy = r * cos(l) * cos(b);
                zz = r * sin(b);
                stars[star_read_num * 3 + 0] = xx;
                stars[star_read_num * 3 + 1] = yy;
                stars[star_read_num * 3 + 2] = zz;
                records[star_read_num * 3 + 0] = l;
                records[star_read_num * 3 + 1] = b;
                records[star_read_num * 3 + 2] = r;
                star_read_num++;
                min_x = xx < min_x ? xx : min_x;
                max_x = xx > max_x ? xx : max_x;
                min_y = yy < min_y ? yy : min_y;
                max_y = yy > max_y ? yy : max_y;
                min_z = zz < min_z ? zz : min_z;
                max_z = zz > max_z ? zz : max_z;
            }
            fclose(fp);
        }
    }
    //print logs
    printf("Stars succesfully loaded.\n");
    printf("\tstars exist: %d\n", star_num);
    printf("\tstars read: %d\n", star_read_num);
    printf("\tmin / max x val: %lf / %lf\n", min_x, max_x);
    printf("\tmin / max y val: %lf / %lf\n", min_y, max_y);
    printf("\tmin / max z val: %lf / %lf\n", min_z, max_z);

    //step-3 generate random (x y z) for kmeans algorithm
    for(int cnt = 0; cnt < centroid_num; ++cnt)
    {
        srand((unsigned int)stars[(2*cnt) % star_num]);
        centroid_x[cnt] = ((rand() % 100) / 100.0)*(max_x - min_x) + min_x;
        srand((unsigned int)stars[(20*cnt) % star_num]);
        centroid_y[cnt] = ((rand() % 100) / 100.0)*(max_y - min_y) + min_y;
        srand((unsigned int)stars[(200*cnt) % star_num]);
        centroid_z[cnt] = ((rand() % 100) / 100.0)*(max_z - min_z) + min_z;
    }
    printf("K-Means Centroid Init:\n");
    for(int cnt = 0; cnt < centroid_num; ++cnt)
    {
        printf("\tc%d(%.5lf,%.5lf,%.5lf) ", cnt, centroid_x[cnt], centroid_y[cnt], centroid_z[cnt]);
    }

    //step-4 call non-parallel k-means clustering
    lloyd_kmeans(stars, star_num, centroid_num, centroid_x, centroid_y, centroid_z, cluster);

    //write non-mpi kmeans result to files
    mkdir("./clusters_non_mpi", 0777);
    FILE ** fa = (FILE **)malloc(sizeof(FILE*)*centroid_num);
    char fn[128] = {'\0'};
    //create and open files
    for(int cnt = 0; cnt < centroid_num; ++cnt)
    {
        sprintf(fn, "./clusters_non_mpi/stars-%d.txt", cnt);
        fa[cnt] = fopen(fn, "w");
        fprintf(fa[cnt], "%d\n", cluster_count[cnt]);
    }
    //write data by cluster id
    for(int idx = 0; idx < star_num; ++idx)
    {
        fprintf(fa[cluster[idx]], "%f %f %f\n", records[idx*3+0], records[idx*3+1], records[idx*3+2]);
    }
    //close files
    for(int cnt = 0; cnt < centroid_num; cnt++)
    {
        fclose(fa[cnt]);
    }
    free(fa);

    //dellocate memory
    free(stars); stars = NULL;
    free(records); records = NULL;
    free(centroid_x); centroid_x = NULL;
    free(centroid_y); centroid_y = NULL;
    free(centroid_z); centroid_z = NULL;
    free(centroid_x_itr); centroid_x_itr = NULL;
    free(centroid_y_itr); centroid_y_itr = NULL;
    free(centroid_z_itr); centroid_z_itr = NULL;
    free(centroid_n_itr); centroid_n_itr = NULL;
    free(cluster); cluster = NULL;
    free(cluster_count); cluster_count = NULL;

    return 0;
}




