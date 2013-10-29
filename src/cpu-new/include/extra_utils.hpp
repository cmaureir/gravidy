#ifndef EXTRAUTILS_HPP
#define EXTRAUTILS_HPP
#include <cmath>
#include "common.hpp"

float  get_magnitude(float x, float y, float z);
double get_timestep_normal(int i, float4 *a2, float4 *a3, double *dt, Forces *f, float eta);
double normalize_dt(double new_dt, double old_dt, double t, int i);

#endif
