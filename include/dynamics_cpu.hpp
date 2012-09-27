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
float initial_energy();
void get_energy_log(int,float,int, float*);
void save_old();
void predicted_pos_vel(float);
void predicted_pos_vel_kepler(float, int);
void correction_pos_vel(float,int);
