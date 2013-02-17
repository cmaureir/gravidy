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
double int_time;
double ini_time, end_time;
double init_time;
double energy_ini, energy_end, energy_tmp;
double ekin, epot;
float softening, eta;

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
/*
 * Main
 */
int
main(int argc, char *argv[])
{
    // Read parameters
    if(!check_options(argc,argv)) return 1;

    read_input_file(input_file);
    alloc_vectors_cpu();

    // Opening output file for debugging
    out = fopen(output_file.c_str(), "w");

    ini_time = (float)clock()/CLOCKS_PER_SEC;
    integrate_cpu();
    end_time = (float)clock()/CLOCKS_PER_SEC;

    fclose(out);
    write_output_file(output_file);
    free_vectors_cpu();

    return 0;
}
