#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <ctime>
#include <omp.h>
#include <cassert>
//#include <thrust/host_vector.h>

/*********************************
 *  Preprocessor directives
 *********************************/

/*
 * Debugging options
 */
//#define KERNEL_ERROR_DEBUG 1
//#define DEBUG_KEPLER 1
//#define DEBUG_HERMITE 1

/*
 * Variable that define is we use the Keplerian correction or not.
 */
//#define USE_KEPLER 1

/*
 * Keplerian correction variables
 */
#ifdef USE_KEPLER
    #define INIT_PARTICLE 1  // Starting from 1 to avoid including the BH
                             // in all the loops
#else
    #define INIT_PARTICLE 0  // Starting from 0 to include all the particles
                                // in all the loops
#endif

#define KEPLER_ITE (50)      // Maximum of iteration when solving Kepler's equation
#define DEL_E      (9.0e-16) // Maximum error in E in elliptical orbits.
#define DEL_E_HYP  (2.e-15)  // Maximum error in E in hyperbolic orbits.
#define OSTEPS     (50)      // Maximum of steps when calculating the central
                             // time-steps distribution

#define J 10 // Neighbour amount to calculate the center of density

/*
 * Softening parameter
 * (Please note that this parameter can be modified by the command line)
 */
#define E 1e-4
#define E2 1e-8

/*
 * ETA_N used to obtain the new time-step for a particle
 * using the equation (7) from Makino and Aarseth 1992.
 * (Please note that this parameter can be modified by the command line)
 *
 * ETA_S used to obtain the firsts time-steps for all the
 * particles of the system.
 */
#define ETA_S 0.001
#define ETA_N 0.01

/*
 * Time-step limits, used to restrict the values
 * using some minimum and maximum boundaries.
 *
 * The limits are powers of 2, because we use the block time-steps scheme.
 */
#define D_TIME_MIN (1.1920928955078125e-07) // 2e-23
#define D_TIME_MAX (0.125)                  // 2e-3

/*
 * Gravitational constant in N-body units
 */
#define G 1

/*
 * CUDA Configuration
 */
#define BSIZE   64  // Block size on kernels calls
#define NJBLOCK 16  // Block size of the shared memory loading j-particles

// Macro from cutil.h to debug the CUDA calls
#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                          \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
    } } while (0)
#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);

// Defining the data type double4 when working on systems without GPU
//typedef struct double4
//{
//    double x, y, z, w;
//} double4;

/*****************************************
 * General variables of the integrator
 *****************************************/

/*
 *  Particle structure to read the input file
 *  (double4 is used instead of double3 to copy the data easily to the
 *  real integrator arrays)
 */
typedef struct particle
{
    float m;
    double4 r;
    double4 v;
} particle;

extern std::vector<particle> part; // Vector to save the input file data


/*
 * General variables of the program
 */
extern int n;                         // Number of particles on the system
extern int iterations;                // Number of iterations in the integration
extern std::string input_file;        // Input filename
extern std::string output_file;       // Output filename
extern std::string run;               // Run/Execution option (cpu or gpu)
extern float total_mass;              // Total mass of the particles
                                      // (In N-body units will be 1)
extern double int_time;               // Integration clock time
extern double ini_time, end_time;     // Initial and Final clock time stamps
extern double ekin, epot;             // Kinetic and Potential energy
extern double energy_ini;             // Initial energy of the system
extern double energy_end;             // Energy at an integration time t
extern double energy_tmp;             // Energy at an integration time t-1
extern float  softening, eta;         // Softening and ETA parameters
                                      // (This parameters takes the previous
                                      // setted values E, and ETA_N or the
                                      // parameters gived by the command line)
extern float t_rh;                    // Half-mass relaxation time
extern float t_cr;                    // Crossing time
extern size_t d1_size, d4_size;       // double and double4 size
extern size_t f1_size, i1_size;       // float and int size
extern size_t nthreads, nblocks;      // Dynamical number of threads and blocks
                                      // for the GPU use.

/*********************************
 *  Host and Devive pointers
 *********************************/

/*
 * Host pointers
 * (Particles attribute arrays)
 */
extern double *h_ekin, *h_epot; // Kinetic and Potential energy
extern double  *h_t, *h_dt;     // Time and timestep
extern double4 *h_r, *h_v;      // Position and Velocity
extern double4 *h_a, *h_a1;     // Acceleration and its first derivative (Jerk)
extern float   *h_m;            // Masses of the particles
extern int *h_move;             // Particles id to move in each iteration time
extern double4 *h_a2, *h_a3;    // 2nd and 3rd acceleration derivatives.
extern double4 *h_old_a;        // Previous step value of the Acceleration
extern double4 *h_old_a1;       // Previous step value of the Jerk
extern double4 *h_p_r;          // Predicted Position
extern double4 *h_p_v;          // Predicted Velocity

/*
 * Device pointers
 * (Particles attribute arrays)
 */
extern double *d_ekin, *d_epot; // Kinetic and Potential energy
extern double  *d_t, *d_dt;     // Time and timestep
extern double4 *d_r, *d_v;      // Position and Velocity
extern double4 *d_a, *d_a1;     // Acceleration and its first derivative (Jerk)
extern float   *d_m;            // Masses of the particles
extern int *d_move;             // Particles id to move in each iteration time
extern double4 *d_p_r;          // Predicted Position
extern double4 *d_p_v;          // Predicted Velocity



/************************************************
 * Special operators for the double4 data type
 ************************************************/
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

inline __host__ __device__ void operator+=(volatile double4 &a, volatile double4 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;

}
