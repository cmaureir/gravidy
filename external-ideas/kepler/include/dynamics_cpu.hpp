#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include "extra_utils.hpp"
#include "kepler.hpp"

void next_itime(double*);
int find_particles_to_move(double);
void update_iteration_dt();
void init_dt(double*);
void init_acc_jrk();
void update_acc_jrk(int);
double energy();
void save_old(int);
void predicted_pos_vel(double);
void correction_pos_vel(double,int);
double get_timestep_normal(int);
double get_timestep_central(int);
double normalize_dt(double, double, double, int);
