#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <thrust/host_vector.h>
#include <omp.h>

// Initial time step factor
#define ETA_S 0.01
// Update time step factor
#define ETA_N 0.02
// Softening parameter
#define E 1e-4
// Block time steps minimum and maximum
#define D_TIME_MIN (10e-8)
#define D_TIME_MAX (0.125)

// Gravitational constant in nbu
#define G 1

// CUDA Fixed block size
#define BSIZE 64

// Neighbor number
#define J 10

#define PI 3.141592653589793

#define MIN_PART_TO_MOVE 10


//typedef struct float4
//{
//    float x,y,z,w;
//} float4;
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
extern float int_time;               // Integration time
extern float ini_time, end_time;     // Initial and Final time stamps
extern float init_time;              // Initial Acc and Jerk calculation time.
extern float energy_ini, energy_end, energy_total; // Initial and Final energy of the system
extern float ekin, epot;             // Kinetic and Potential energy

/*
 * Struct vector to read the input file
 */
extern std::vector<particle> part;


/*
 * Host pointers
 * Particles attribute arrays
 */
extern float *h_ekin, *h_epot;         // Kinetic and Potential energy
extern float *h_t, *h_dt;              // Time and timestep
extern double4 *h_r, *h_v, *h_a, *h_j;  // Position, velocity, acceleration and jerk
extern float   *h_m;                    // Masses
extern int *h_move;                     // Index of the particles to move in each iteration time
extern double4 *h_new_a, *h_new_j;      // Temp arrays to save tmp accelerations

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
extern float *d_ekin, *d_epot;         // Kinetic and Potential energy
extern float *d_t, *d_dt;              // Time and timestep
extern double4 *d_r, *d_v, *d_a, *d_j;  // Position, velocity, acceleration and jerk
extern float   *d_m;                    // Masses
extern int *d_move;                     // Index of the particles to move in each iteration time
extern double4 *d_new_a, *d_new_j;      // Temp arrays to save tmp accelerations


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
