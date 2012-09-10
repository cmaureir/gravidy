#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

#include <vector>
#include <cmath>
#include <algorithm>
typedef struct distance
{
    int index;
    double value;
    bool operator<(const struct distance& a) const
    {
        return value < a.value;
    }
} distance;

typedef struct point
{
    double x,y,z;
} point;

point get_center_of_density();
double get_halfmass_radius(double, double, double);

float get_crossing_time();
float get_relaxation_time();
