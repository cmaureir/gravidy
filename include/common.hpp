#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <ctime>
//#include <thrust/host_vector.h>
#include <omp.h>
#include <cassert>

//#define KERNEL_ERROR_DEBUG 1
//#define USE_KEPLER 1
//#define DEBUG_KEPLER 1
//#define DEBUG_HERMITE 1

// The following values can be modified by the command
// line parameters
// Initial time step factor
#define ETA_S 0.001
// Update time step factor
#define ETA_N 0.01 // Using 0.1 to work with a BH
// Softening parameter
#define E 1e-4
#define E2 1e-8
// Block time steps minimum and maximum
#define D_TIME_MIN (1.1920928955078125e-07) // 2e-23
#define D_TIME_MAX (0.125) // 2e-3

// NBODY6 limits
//#define D_TIME_MIN (3.0517578125e-05)
//#define D_TIME_MAX (1)

// Gravitational constant in nbu
#define G 1

// CUDA Fixed block size
#define BSIZE 64
#define NJBLOCK 16


#ifdef USE_KEPLER
#define INIT_PARTICLE    1
#else
#define INIT_PARTICLE    0
#endif

// Neighbor number
#define J                10
#define KEPLER_ITE       50
#define DEL_E            9.0e-16 // Max error in E kepler equation
#define DEL_E_HYP        2.e-15  // Max error in E for hyperbolic kepler equation
#define NSTEPS           50


// Macro from cutil.h
#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                          \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
    } } while (0)
#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);

//typedef struct double4
//{
//    double x, y, z, w;
//} double4;

// Particle structure to read input file
typedef struct particle
{
    float m;
    double4 r;
    double4 v;
} particle;

/*
 * General variables of the program
 */
extern int n;                         // Number of particles
extern int iterations;                // Number of iterations
extern float total_mass;              // Total mass of the particles
extern double int_time;               // Integration time
extern double ini_time, end_time;     // Initial and Final time stamps
extern double init_time;              // Initial Acc and Jerk calculation time.
extern double energy_ini, energy_end, energy_tmp; // Initial and Final energy of the system
extern double ekin, epot;             // Kinetic and Potential energy
extern float softening, eta;

// Struct vector to read the input file
extern std::vector<particle> part;


/*
 * Host pointers
 * Particles attribute arrays
 */
extern double *h_ekin, *h_epot;         // Kinetic and Potential energy
extern double  *h_t, *h_dt;              // Time and timestep
extern double4 *h_r, *h_v, *h_a, *h_a1;  // Position, velocity, acceleration and jerk
extern float   *h_m;                    // Masses
extern int *h_move;                     // Index of the particles to move in each iteration time
//extern double4 *h_new_a, *h_new_j;      // Temp arrays to save tmp accelerations
extern double4 *h_a2, *h_a3;

/*
 * Host pointers
 * Previous step attributes
 */
extern double4 *h_old_a; // Acceleration
extern double4 *h_old_a1; // Jerk

/*
 * Host pointers
 * Predicted attributes
 */
extern double4 *h_p_r; // Position
extern double4 *h_p_v; // Velocity

/*
 * Device pointers
 * Particles attribute arrays
 */
extern double *d_ekin, *d_epot;         // Kinetic and Potential energy
extern double *d_t, *d_dt;              // Time and timestep
extern double4 *d_r, *d_v, *d_a, *d_a1; // Position, velocity, acceleration and jerk
extern float   *d_m;                   // Masses
extern int *d_move;                    // Index of the particles to move in each iteration time
//extern double4 *d_new_a, *d_new_j;     // Temp arrays to save tmp accelerations


/*
 * Device pointers
 * Previous step attributes
 */
extern double4 *d_old_a; // Acceleration
extern double4 *d_old_a1; // Jerk

/*
 * Device pointers
 * Predicted attributes
 */
extern double4 *d_p_r; // Position
extern double4 *d_p_v; // Velocity

/*
 * General system time
 */
extern float t_rh; // Half-mass relaxation time
extern float t_cr; // Crossing time

/*
 * Options strings
 */
extern std::string input_file, output_file; // Input and Output filename
extern std::string run;                     // Run option (cpu or gpu)

extern size_t nthreads, nblocks;
extern size_t d1_size, d4_size;
extern size_t f1_size, i1_size;
/*
 * General system time¬
 */
extern float t_rh; // Half-mass relaxation time¬
extern float t_cr; // Crossing time¬

extern double4 *tmp_red;

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
