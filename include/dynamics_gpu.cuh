#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include "extra_utils.hpp"
#include "dynamics_gpu_kernels.cuh"

__host__ void   gpu_init_dt(float*);
__host__ void   gpu_next_itime(float*);
__host__ int    gpu_find_particles_to_move(float);
__host__ double gpu_energy();
__host__ void   gpu_predicted_pos_vel(float);
__host__ void   gpu_init_acc_jrk();
__host__ void   gpu_update_acc_jrk(int);
__host__ void   gpu_update_2d(int);
__host__ void   gpu_update_acc_jrk_single(int);
__host__ void   gpu_correction_pos_vel(float, int);
