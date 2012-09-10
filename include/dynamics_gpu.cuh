#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include "dynamics_gpu_kernels.cuh"

__host__ double gpu_energy(bool);
__host__ void   gpu_init_acc_jerk();
__host__ void   gpu_update_acc_jerk_simple(int);
__host__ void   gpu_update_acc_jerk_tile(int);
__host__ void   gpu_update_acc_jerk_single(int);
__host__ void   gpu_correction_pos_vel(double, int);
