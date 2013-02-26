#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include "extra_utils.hpp"
#include "dynamics_gpu_kernels.cuh"

__host__ void gpu_get_energy_log(double,int, int);

__host__ double gpu_energy();
__host__ void   gpu_init_acc_jrk();
__host__ void gpu_update_acc_jrk_simple(int);

__host__ void gpu_update(int);
