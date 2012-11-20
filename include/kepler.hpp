#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif
#include <cmath>
#include <iostream>

void kepler_prediction(double*, double*, double*,
                       double*, double*, double*,
                       double, int);
double solve_kepler(double, double);
double kepler(const double, double);
