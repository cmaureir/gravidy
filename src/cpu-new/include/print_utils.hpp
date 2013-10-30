#ifndef PRINTUTILS_HPP
#define PRINTUTILS_HPP
#include "common.hpp"
#include <iostream>
#include <cstdio>
#include <cmath>


void print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f, double *dt);
void print_energy_log(double ITIME, int iterations, long long interactions,
                      int nsteps, Gtime gtime, Energy &energy, double new_energy);
#endif
