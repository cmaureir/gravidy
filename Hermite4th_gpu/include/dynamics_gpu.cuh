#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include "extra_utils.hpp"
#include "dynamics_gpu_kernels.cuh"

__host__ void gpu_get_energy_log(double,int, int);
__host__ void gpu_send_initial_data();

__host__ void gpu_recv_initial_data();
__host__ double gpu_energy();
__host__ void   gpu_predicted_pos_vel(float);
__host__ void   gpu_init_acc_jrk();
__host__ void   gpu_update_acc_jrk(int);
__host__ void   gpu_update_2d(int);
__host__ void   gpu_update_acc_jrk_single(int);
__host__ void gpu_correction_pos_vel(double,int);
__host__ void gpu_update_acc_jrk_simple(int);
__host__ double normalize_dt(double, double, double, int);
