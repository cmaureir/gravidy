#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include "extra_utils.hpp"

void next_itime(double*);
int find_particles_to_move(double);
void init_dt(double*);
void save_old(int);
double get_timestep_normal(int);
double normalize_dt(double, double, double, int);
