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
void print_all(int, float);
void print_positions(int);
void print_velocities(int);
void print_accelerations(int);
void print_jrks(int);
void print_accjrk(int);
void print_masses(int);
void print_times(int, float);
void print_old(int);
void print_predicted(int);
void print_movement(int, int, float);

// Single
void print_particle(int);

// Functions
double get_magnitude(double, double, double);
void get_energy_log(double,int, int, FILE*, double);
