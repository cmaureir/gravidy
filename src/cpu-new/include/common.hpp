#ifndef COMMON_HPP
#define COMMON_HPP
#include <iomanip>
#include <omp.h>

#include <iostream>

const int J=10;
const int INIT_PARTICLE=0;
const float RADIUS_MASS_PORCENTAGE=0.2;
const double E=1e-4;
const double E2=1e-8;
const double  ETA_S=0.01;
const double  ETA_N=0.01;
const double D_TIME_MIN=1.1920928955078125e-07;
const double D_TIME_MAX=0.125;

typedef struct Energy
{
    double ini;
    double end;
    double tmp;
    double kinetic;
    double potential;
} Energy;

typedef struct double4
{
    double x, y, z, w;
} double4;

typedef struct Predictor {
    double r[3];
    double v[3];
    float  m;
} Predictor;

typedef struct Forces {
    double a[3];
    double a1[3];
} Forces;

typedef struct Gtime {
    double integration_ini;
    double integration_end;

    double prediction_ini;
    double prediction_end;

    double update_ini;
    double update_end;

    double correction_ini;
    double correction_end;

    double grav_ini;
    double grav_end;

    double reduce_ini;
    double reduce_end;

    float gflops;
} Gtime;


typedef struct file_data
{
    float m;
    double r[3];
    double v[3];
} file_data;



#endif
