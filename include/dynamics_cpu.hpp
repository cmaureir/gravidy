#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include "extra_utils.hpp"
#include "kepler.hpp"

void next_itime(float*);
int find_particles_to_move(float);
void update_iteration_dt();
void init_dt(float*);
void init_acc_jrk();
void update_acc_jrk(int);
float energy();
void get_energy_log(float,int, int);
void save_old(int);
void predicted_pos_vel(float);
void predicted_pos_vel_kepler(float, int);
void correction_pos_vel(float,int);
float normalize_dt(float, float, float, int);
