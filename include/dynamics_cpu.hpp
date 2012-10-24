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
float energy();
void get_energy_log(double,int, int);
void save_old(int);
void predicted_pos_vel(double);
void predicted_pos_vel_kepler(float, int);
void correction_pos_vel(double,int);
double normalize_dt(double, double, double, int);
