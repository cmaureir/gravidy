#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <thrust/host_vector.h>
#include <omp.h>
#include <cassert>

//#define KERNEL_ERROR_DEBUG 1

// Initial time step factor
#define ETA_S 0.01
// Update time step factor
#define ETA_N 0.01
// Softening parameter
#define E 1e-4
#define E2 1e-8
// Block time steps minimum and maximum
#define D_TIME_MIN (1.1920928955078125e-07) // 2e-23
#define D_TIME_MAX (0.125) // 2e-3

#define ITE_MAX (1e6)
#define OUT (0.01)

// Gravitational constant in nbu
#define G 1

// CUDA Fixed block size
#define BSIZE 64
#define NJBLOCK 16

// Neighbor number
#define J                10
#define PI               3.141592653589793
#define TWOPI            (2*3.14159265)
#define MIN_PART_TO_MOVE 10
#define KEPLER_ITE       50
#define INIT_PARTICLE    0
#define DEL_E            9.0e-16 // Max error in E kepler equation
#define DEL_E_HYP        2.e-15  // Max error in E for hyperbolic kepler equation


// Macro from cutil.h
#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                          \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
    } } while (0)
#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);



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

// Struct vector to read the input file
extern std::vector<particle> part;


/*
 * Host pointers
 * Particles attribute arrays
 */
extern double *h_ekin, *h_epot;         // Kinetic and Potential energy
extern double *h_t, *h_dt;              // Time and timestep
extern double4 *h_r, *h_v, *h_a, *h_j;  // Position, velocity, acceleration and jerk
extern float   *h_m;                    // Masses
extern int *h_move;                     // Index of the particles to move in each iteration time
//extern double4 *h_new_a, *h_new_j;      // Temp arrays to save tmp accelerations

/*
 * Host pointers
 * Previous step attributes
 */
extern double4 *h_old_r; // Position
extern double4 *h_old_v; // Velocity
extern double4 *h_old_a; // Acceleration
extern double4 *h_old_j; // Jerk

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
extern double4 *d_r, *d_v, *d_a, *d_j; // Position, velocity, acceleration and jerk
extern float   *d_m;                   // Masses
extern int *d_move;                    // Index of the particles to move in each iteration time
//extern double4 *d_new_a, *d_new_j;     // Temp arrays to save tmp accelerations


/*
 * Device pointers
 * Previous step attributes
 */
extern double4 *d_old_r; // Position
extern double4 *d_old_v; // Velocity
extern double4 *d_old_a; // Acceleration
extern double4 *d_old_j; // Jerk

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
