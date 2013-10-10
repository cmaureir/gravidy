#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include "extra_utils.hpp"
#include "integration_utils.hpp"

void force_calculation(int, int);
void init_acc_jrk();
void update_acc_jrk(int);
void predicted_pos_vel(double);
void correction_pos_vel(double,int);
double energy();
