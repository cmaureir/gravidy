#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

__global__ void k_energy(double4*,
                         double4*,
                         double*,
                         double*,
                         float*,
                         int);

__global__ void k_init_acc_jrk(double4*,
                               double4*,
                               Forces*,
                               float*,
                               int);

__global__ void k_predicted_pos_vel(double4*,
                                    double4*,
                                    Forces*,
                                    Predictor*,
                                    double*,
                                    double,
                                    int);

__device__ void k_force_calculation(double4,
                                    double4,
                                    double4,
                                    double4,
                                    Forces&,
                                    float);

__device__ void k_force_calculation2(Predictor,
                                     Predictor,
                                     Forces&,
                                     float);

__global__ void k_update(Predictor*,
                         Predictor*,
                         Forces*,
                         float*,
                         int,
                         int);

__global__ void reduce(Forces*,
                       Forces*,
                       unsigned int);
