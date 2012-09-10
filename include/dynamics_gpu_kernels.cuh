#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

__global__ void k_energy                ( double4*, double4*, double*,  double*,  float*, int);
__global__ void k_init_acc_jerk         ( double4*, double4*, double4*, double4*, float*, int);
__global__ void k_init_acc_jerk_tile    ( double4*, double4*, double4*, double4*, float*, int);
__global__ void k_update_acc_jerk_simple( double4*, double4*, double4*, double4*, float*, int*, int, int);
__global__ void k_update_acc_jerk_tile  ( double4*, double4*, double4*, double4*, float*, int*, int, int);
__global__ void k_correction_pos_vel    ( double4*, double4*, double4*, double4*,
                                          double4*, double4*,
                                          double4*, double4*, double*,  double*,  
                                          double,   int*,     int);
__global__ void k_update_acc_jerk_single( double4,  double4,  double4*, double4*, double4*,
                                          double4*, float*,   int,      int);
