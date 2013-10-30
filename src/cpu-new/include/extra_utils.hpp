#ifndef EXTRAUTILS_HPP
#define EXTRAUTILS_HPP
#include <cmath>
#include "common.hpp"

double  get_magnitude(double x, double y, double z);
double get_timestep_normal(int i, double4 *a2, double4 *a3, double *dt, Forces *f, float eta);
double normalize_dt(double new_dt, double old_dt, double t, int i);

#endif
