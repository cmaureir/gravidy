#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

struct Distance
{
    int index;
    double value;
    bool operator<(const Distance& a) const
    {
        return value < a.value;
    }
};

struct Point
{
    double x,y,z;
};

Point get_center_of_density();
double get_halfmass_radius(double, double, double);
float get_crossing_time();
float get_relaxation_time();
