#ifndef COMMON_HPP
#define COMMON_HPP
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <omp.h>

#ifndef GPU
typedef struct double4
{
    double x, y, z, w;
} double4;

typedef struct double3
{
    double x, y, z;
} double3;
#else
#include <vector_types.h>
#endif

const int    G                      = 1;
const int    J                      = 10;
const int    INIT_PARTICLE          = 0;
const float  RADIUS_MASS_PORCENTAGE = 0.2;
const double E                      = 1e-4;
const double E2                     = 1e-8;
const double ETA_S                  = 0.01;
const double ETA_N                  = 0.01;
const double D_TIME_MIN             = 1.1920928955078125e-07;
const double D_TIME_MAX             = 0.125;

typedef struct Energy
{
    double ini;
    double end;
    double tmp;
    double kinetic;
    double potential;
} Energy;


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

#ifdef GPU
const int BSIZE   = 32;
const int NJBLOCK = 16;

// Macro from cutil.h to debug the CUDA calls
#define CUDA_SAFE_CALL_NO_SYNC( call) do {                          \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
    } } while (0)
#define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);

inline __host__ __device__ double4 operator+(const double4 &a, const double4 &b)
{
    double4 tmp = {a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w};
    return tmp;
}

inline __host__ __device__ void operator+=(double4 &a, double4 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
}

inline __host__ __device__ Forces operator+(Forces &a, Forces &b)
{
    Forces tmp;
    tmp.a[0] = a.a[0] + b.a[0];
    tmp.a[1] = a.a[1] + b.a[1];
    tmp.a[2] = a.a[2] + b.a[2];

    tmp.a1[0] = a.a1[0] + b.a1[0];
    tmp.a1[1] = a.a1[1] + b.a1[1];
    tmp.a1[2] = a.a1[2] + b.a1[2];

    return tmp;
}

inline __host__ __device__ void operator+=(Forces &a, Forces &b)
{
    a.a[0] += b.a[0];
    a.a[1] += b.a[1];
    a.a[2] += b.a[2];

    a.a1[0] += b.a1[0];
    a.a1[1] += b.a1[1];
    a.a1[2] += b.a1[2];
}
#endif

#endif
