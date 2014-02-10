#ifndef COMMON_HPP
#define COMMON_HPP

#ifdef _MPI
#include <mpi.h>
#endif

#include <iomanip>
#include <iostream>
#include <cstdio>
#include <omp.h>


#ifndef GPU
/** Defining the «double4» structure based on the CUDA definition for the
 * CPU version, which not include the CUDA headers */
typedef struct double4
{
    double x, y, z, w;
} double4;

/** Defining the «double3» structure based on the CUDA definition for the
 * CPU version, which not include the CUDA headers */
typedef struct double3
{
    double x, y, z;
} double3;
#else
/** If we are compiling the CUDA version, we add the definition of the
 * vector types and structs from the CUDA library */
#include <vector_types.h>
#endif

/** Gravitational constant. Since we are working in N-body units we set G as one. */
const int    G                      = 1;
/** Amount of neighbours to calculate the center of density of the system */
const int    J                      = 10;
/** Common radius for the core of a globular cluster */
const float  RADIUS_RATIO = 0.05;
/** Softening parameter */
const double E                      = 1e-4;
/** Softening parameter squared */
const double E2                     = 1e-8;
/** Initial ETA parameter to calculate the first timestep of all the particles
 * of the system. Based on Aarseth formula */
const double ETA_S                  = 0.01;
/** Iteration ETA parameter to calculate new timestep of all the active particles
 * of the system, in a certain integration time. Based on Aarseth formula */
const double ETA_N                  = 0.01;
/** For Kepler **/
const double ETA_K                  = 0.1;
/** Lower boundary for the particles timesteps, \f$2^{-23}\f$ */
const double D_TIME_MIN             = 1.1920928955078125e-07;
/** Upper boundary for the particles timesteps, \f$2^{-3}\f$ */
const double D_TIME_MAX             = 0.125;

#ifdef KEPLER
const int FIRST_PARTICLE = 1;
#else
const int FIRST_PARTICLE = 0;
#endif


#define KEPLER_ITE (50)      // Maximum of iteration when solving Kepler's equation
#define DEL_E      (9.0e-16) // Maximum error in E in elliptical orbits.
#define DEL_E_HYP  (2.e-15)  // Maximum error in E in hyperbolic orbits.
#define OSTEPS     (50)      // Maximum of steps when calculating the central
                             // time-steps distribution
struct Distance
{
    int index;
    double value;
    bool operator<(const Distance& a) const
    {
        return value < a.value;
    }
};

/** @struct Energy
 *  @brief This structure contains all the energy variables of the system.
 *  @var Energy::ini
 *  Member 'ini' contains the initial total energy of the system.
 *  @var Energy::ini
 *  Member 'ini' contains the initial total energy of the system.
 *  @var Energy::end
 *  Member 'end' contains the newest total energy of the system in a certain time.
 *  @var Energy::tmp
 *  Member 'tmp' contains the previous total total energy of the system.
 *  @var Energy::kinetic
 *  Member 'kinetic' contains the newest kinetic energy of the system in a certain time.
 *  @var Energy::potential
 *  Member 'potential' contains the newest kinetic energy of the system in a certain time.
 *  */
typedef struct Energy
{
    double ini;
    double end;
    double tmp;
    double kinetic;
    double potential;
} Energy;

typedef struct options
{
    bool print_screen;
    bool print_all;
    bool print_lagrange;
} options;

/** @struct Predictor
 *  @brief This structure contains the predicted information of a particle in some
 *         moment of the integration.
 *  @var Predictor::r
 *  Member 'r' contains the position in three dimensions.
 *  @var Predictor::v
 *  Member 'v' contains the velocity in three dimensions.
 *  @var Predictor::m
 *  Member 'm' contains the mass of the particle.
 *  */
typedef struct Predictor {
    double r[3];
    double v[3];
    float  m;
} Predictor;

/** @struct Forces
 *  @brief This structure contains the information of the Forces of a particle in some
 *         moment of the integration.
 *  @var Forces::a
 *  Member 'a' contains the acceleration in three dimensions.
 *  @var Forces::a1
 *  Member 'v' contains the first derivative of the acceleration in three dimensions
 *  (Jerk).
 *  */
typedef struct Forces {
    double a[3];
    double a1[3];
} Forces;

/** @struct Gtime
 *  @brief This structure contains different times of the internal integration
 *         process.
 *         This times are calculated using the function omp_get_wtime() from the
 *         OpenMP library.
 *  @var Gtime::integration_ini
 *  Member 'integration_ini' contains the starting time of the integration.
 *  @var Gtime::integration_end
 *  Member 'integration_end' contains the final time of the integration.
 *  @var Gtime::prediction_ini
 *  Member 'prediction_ini' contains the starting time of the prediction.
 *  @var Gtime::prediction_end
 *  Member 'prediction_end' contains the final time of the prediction.
 *  @var Gtime::update_ini
 *  Member 'update_ini' contains the starting time of the forces update.
 *  @var Gtime::update_end
 *  Member 'update_end' contains the final time of the forces update.
 *  @var Gtime::correction_ini
 *  Member 'correction_ini' contains the starting time of the correction.
 *  @var Gtime::correction_end
 *  Member 'correction_end' contains the final time of the correction.
 *  @var Gtime::grav_ini
 *  Member 'grav_ini' contains the starting time of the gravitational interaction.
 *  @var Gtime::grav_end
 *  Member 'grav_end' contains the final time of the gravitational interaction.
 *  @var Gtime::reduce_ini
 *  Member 'reduce_ini' contains the starting time of the forces reduction.
 *  @var Gtime::reduce_end
 *  Member 'reduce_end' contains the final time of the forces reduction.
 *  @var Gtime::gflops
 *  Member 'gflops' contains the amount of Giga FLOPs of the force update method.
    This is calculated with the following formula:
        \f$ 60.10e-9 \cdot \frac{1}{C_{\rm time}}\cdot \sum_{t=0}^{t=T} N_{\rm act} N \f$
    where \f$(N_{\rm act} N)\f$ is the amount of gravitational interactions,
    \f$C_{\rm time}\f$ the elapsed clock-time of the process,
    \f$T\f$ a determinated integration time.
 *  */
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
