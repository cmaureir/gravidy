#ifndef COMMON_HPP
#define COMMON_HPP
#include <iostream>
#include <iomanip>

const int J=10;
const int INIT_PARTICLE=0;
const float RADIUS_MASS_PORCENTAGE=0.2;
const float E=1e-4;
const float E2=1e-8;
const float ETA_S=0.01;
const float ETA_N=0.01;
const float D_TIME_MIN=1.1920928955078125e-07;
const float D_TIME_MAX=0.125;

typedef struct double4
{
    double x, y, z, w;
} double4;

typedef struct float4
{
    float x, y, z, w;
} float4;

typedef struct Predictor {
    double r[3];
    double v[3];
} Predictor;

typedef struct Forces {
    float a[3];
    float a1[3];
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
} Gtime;


typedef struct file_data
{
    float m;
    float r[3];
    float v[3];
} file_data;

#endif
