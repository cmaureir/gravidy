#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

__global__ void k_update_2d             ( int*, double4*, double4*, double4*, double4*, float*, int, int);
__global__ void k_reduce_energy         ( float*, float*, int);
__global__ void k_reduce                ( double4*, double4*, int);
__global__ void k_energy                ( double4*, double4*, float*,   float*,   float*, int);
__global__ void k_init_acc_jrk         ( double4*, double4*, double4*, double4*, float*, int);
__global__ void k_update_acc_jrk       ( double4*, double4*, double4*, double4*, float*, int*, int, int);
__global__ void k_correction_pos_vel    ( double4*, double4*, double4*, double4*,
                                          double4*, double4*,
                                          double4*, double4*, float*,  float*,  
                                          float,    int*,     int);
__global__ void k_update_acc_jrk_single( double4,  double4, double4*, double4*, double4*,
                                          double4*, float*,   int,      int);

__global__ void
k_predicted_pos_vel(double4*, double4*, double4*, double4*,
                    double4*, double4*, float*,   float, int);

__global__ void
k_update_acc_jrk_simple(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int *move, int n, int total);
