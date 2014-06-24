#ifndef NBODYSYSTEM_HPP
#define NBODYSYSTEM_HPP

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <algorithm>

#include "common.hpp"
#include "utils/OptionsParser.hpp"

/*
 * @class NbodySystem
 * @desc N-body system general class with all the information of the integration
 *       process and the particles of the system.
 */
class NbodySystem {
    public:
        NbodySystem(OptionsParser op);
        ~NbodySystem();

        Gtime gtime;
        std::vector<file_data> reader;
        options ops;

        // Configuration parameters
        int    n;
        float  eta;
        float  total_mass;
        float  integration_time;
        long long int iterations;
        double e2;

        // Global parameters
        Energy en;
        float  r_h;
        float  hmr_time;
        float  hmr_time_soft;
        float  cr_time;
        double3 cod;

        // Files
        std::string   input_filename;
        std::string   output_filename;

        /******************************** Host Particles arrays */

        int       *h_move; // Array with index of the active particles
        double    *h_t;    // Particle's time
        double    *h_dt;   // Particle's timestep
        double    *h_ekin; // Kinetic energy
        double    *h_epot; // Potential Energy
        double4   *h_r;    // Position + Mass
        double4   *h_v;    // Velocity + Empty
        double3   *h_a2;   // 2nd acceleration derivative (from interpolation)
        double3   *h_a3;   // 3rd acceleration derivative (from interpolation)
        Predictor *h_p;    // Predicted position and velocities
        Predictor *h_i;    // Temporary array of the position and velocity of the
                           //   active particle only
        Forces    *h_f;    // Struct with the acceleration and jerk
        Forces    *h_old;  // Old forces (from previous step)
        Forces    *h_fout_tmp; // Temporary array for partial forces

        /******************************** Device Particles arrays */

        int       *d_move;
        double    *d_t;
        double    *d_dt;
        double    *d_ekin;
        double    *d_epot;
        double4   *d_r;
        double4   *d_v;
        Predictor *d_p;
        Predictor *d_i;
        Forces    *d_f;
        Forces    *d_fout;
        Forces    *d_fout_tmp;
        Forces    *d_old;

        /******************************** Keplerian correction */
        double  *h_tlast;
        double  *h_dtlast;
        double3 *h_a2bh;
        double3 *h_a3bh;
        Forces  *h_fbh;

        /******************************** General functions of the system */
        void read_input_file();
        void copy_input_data();
        void alloc_base_attributes(int rank);
        void free_base_attributes();
};
#endif // NBODYSYSTEM_HPP
