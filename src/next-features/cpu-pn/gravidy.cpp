#include "include/options_parser.hpp"
#include "include/file_utils.hpp"
#include "include/memory_cpu.hpp"
#include "include/hermite.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <omp.h>

// General variables of the program
int n;
int iterations;
float total_mass;
double energy_ini, energy_end, energy_tmp;
double ekin, epot;
float e2, eta;
Gtime gtime;
float itime;
float gflops;

// Struct vector to read the input file
std::vector<particle> part;
double4  *h_old_a, *h_old_a1;
Predictor *h_p;

// Host pointers
double *h_ekin, *h_epot;
double *h_t, *h_dt;
float   *h_m;
int *h_move;
double4 *h_r, *h_v;
Forces *h_f;
double4 *h_a2, *h_a3;

// System times
float t_rh, t_cr;

// Options strings
std::string input_file, output_file;
FILE *out;

size_t d1_size, d4_size;
size_t f1_size, i1_size;

int print_log;

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

    // Opening output file for debugging
    if (print_log)
    {
        output_file += "_";
        output_file += get_time();
        output_file += ".out.gpu";
        out = fopen(output_file.c_str(), "w");
    }


    // Start integration process
    gtime.integration_ini = omp_get_wtime();
    integrate_cpu();
    gtime.integration_end = omp_get_wtime();

    if(print_log)
    {
        fclose(out);
    }

    // Memory deallocation of the GPU arrays
    free_vectors_cpu();

    return 0;
}
