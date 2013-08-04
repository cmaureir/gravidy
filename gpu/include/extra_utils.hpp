#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

// Group of particles
void print_all(int, float, FILE*);
void print_forces(int);
void print_times(int, float);
void print_predicted(int);
void print_movement(int, int, float);

// Single
void print_particle(int);

// Functions
double get_magnitude(double, double, double);
void get_energy_log(double,int, int, FILE*, double);
string get_time();
void get_kernel_error();
void gpu_timer_start();
float gpu_timer_stop(string);
