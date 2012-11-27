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
#define RADIUS_MASS_PORCENTAGE 0.4

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

typedef struct double4
{
    double x, y, z, w;
} double4;

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
extern std::string output_file;       // Output filename for general info.
extern FILE *out;                     // Out file for debugging.
extern std::string run;               // Run/Execution option (CPU or GPU)
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
                                      // parameters give by the command line)
extern float t_rh;                    // Half-mass relaxation time
extern float t_cr;                    // Crossing time
extern size_t d1_size, d4_size;       // double and double4 size
extern size_t f1_size, i1_size;       // float and int size

/*********************************
 *  Host and Device pointers
 *********************************/

/*
 * Host pointers
 * (Particles attribute arrays)
 */
extern double *h_ekin, *h_epot; // Kinetic and Potential energy
extern double  *h_t, *h_dt;     // Time and time-step
extern double4 *h_r, *h_v;      // Position and Velocity
extern double4 *h_a, *h_a1;     // Acceleration and its first derivative (Jerk)
extern float   *h_m;            // Masses of the particles
extern int *h_move;             // Particles id to move in each iteration time
extern double4 *h_a2, *h_a3;    // 2nd and 3rd acceleration derivatives.
extern double4 *h_old_a;        // Previous step value of the Acceleration
extern double4 *h_old_a1;       // Previous step value of the Jerk
extern double4 *h_p_r;          // Predicted Position
extern double4 *h_p_v;          // Predicted Velocity
