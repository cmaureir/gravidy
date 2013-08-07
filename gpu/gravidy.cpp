#include "include/options_parser.hpp"
#include "include/file_utils.hpp"
#include "include/memory_cpu.hpp"
#include "include/memory_gpu.cuh"
#include "include/hermite.cuh"
#include <time.h>

#include <iostream>
#include <iomanip>
#include <limits>

// General variables of the program
int n;
int iterations;
float total_mass;
double energy_ini, energy_end, energy_tmp;
double ekin, epot;
float e2, eta;
float beta;
Gtime gtime;
float itime;
float gflops;

// Struct vector to read the input file
std::vector<particle> part;
double4  *h_old_a, *h_old_a1;
Predictor *h_p;
Predictor *d_p;

// Host pointers
double *h_ekin, *h_epot;
double *h_t, *h_dt;
float   *h_m;
int *h_move;
double4 *h_r, *h_v;
Forces *h_f;
double4 *h_a2, *h_a3;
Predictor *h_i;
Forces *h_fout;
Forces *h_fout_tmp;

// Device pointers
double *d_ekin, *d_epot;
double  *d_t, *d_dt;
double4 *d_r, *d_v;
Forces *d_f;
float   *d_m;
int *d_move;
Predictor *d_i;
Forces *d_fout;
Forces *d_fout_tmp;

// System times
float t_rh, t_cr;

// Options strings
std::string input_file, output_file;
FILE *out;

size_t d1_size, d4_size;
size_t f1_size, i1_size;
size_t nthreads, nblocks;

int print_log;
cudaEvent_t start, stop;

/********************************************************
 * Main
 ********************************************************/

int
main(int argc, char *argv[])
{
    // Read command line parameters
    if(!check_options(argc,argv)) return 1;

    // Read the input file
    read_input_file(input_file);

    // Memory allocation of the CPU arrays
    alloc_vectors_cpu();

    // Memory allocation of the GPU arrays
    alloc_vectors_gpu();

    // Opening output file for debugging
    if (print_log)
    {
        output_file += "_";
        output_file += get_time();
        output_file += ".out.gpu";
        out = fopen(output_file.c_str(), "w");
    }

    // Create GPU timers
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Start integration process
    gtime.integration_ini = omp_get_wtime();
    integrate_gpu();
    gtime.integration_end = omp_get_wtime();

    if(print_log)
    {
        fclose(out);
    }

    // Memory deallocation of the CPU arrays
    free_vectors_gpu();

    // Memory deallocation of the GPU arrays
    free_vectors_cpu();

    return 0;
}
