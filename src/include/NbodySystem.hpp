/*
 * Copyright (c) 2016
 *
 * Cristián Maureira-Fredes <cmaureirafredes@gmail.com>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef NBODYSYSTEM_HPP
#define NBODYSYSTEM_HPP

#include <iterator>
#include <algorithm>

#include "common.hpp"
#include "utils/OptionsParser.hpp"
#define MAXGPUS 4

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
        int  gpus;
        long long int iterations;
        unsigned int n;
        unsigned int snapshot_number;

        float  eta;
        float  total_mass;
        float  integration_time;

        double e2;
        double snapshot_time;
        bool resume;

        // General variables
        Energy en;
        float t_rlx;
        float t_cr;
        float r_virial;
        float max_mass;

        // Close encounter
        double r_cl;
        double dt_cl;
        double r_core;

        // Heaviest star
        double m_g;

        // Files
        std::string   input_filename;
        std::string   output_filename;
        std::string   resume_filename;
        std::string   snapshot_filename;

        /******************************** Host Particles arrays */

        unsigned int  *h_id;   // Array with index of the particles
        unsigned int  *h_move; // Array with index of the active particles
        double        *h_t;    // Particle's time
        double        *h_dt;   // Particle's timestep
        double        *h_ekin; // Kinetic energy
        double        *h_epot; // Potential Energy
        double4       *h_r;    // Position + Mass
        double4       *h_v;    // Velocity + Empty
        double3       *h_a2;   // 2nd acceleration derivative (from interpolation)
        double3       *h_a3;   // 3rd acceleration derivative (from interpolation)
        Predictor     *h_p;    // Predicted position and velocities
        Predictor     *h_i;    // Temporary array of the position and velocity of the
                               //   active particle only

        Forces        *h_f;                 // Struct with the acceleration and jerk
        Forces        *h_old;               // Old forces (from previous step)
        Forces        *h_fout_tmp;          // Temporary array for partial forces
        Forces        *h_fout_gpu[MAXGPUS]; // Temporary array for partial forces

        /******************************** Device Particles arrays */

        unsigned int *d_id[MAXGPUS];
        unsigned int *d_move[MAXGPUS];
        double       *d_t[MAXGPUS];
        double       *d_dt[MAXGPUS];
        double       *d_ekin[MAXGPUS];
        double       *d_epot[MAXGPUS];
        double4      *d_r[MAXGPUS];
        double4      *d_v[MAXGPUS];
        Predictor    *d_p[MAXGPUS];
        Predictor    *d_i[MAXGPUS];
        Forces       *d_f[MAXGPUS];
        Forces       *d_fout[MAXGPUS];
        Forces       *d_fout_tmp[MAXGPUS];
        Forces       *d_old[MAXGPUS];

        /******************************** General functions of the system */
        void read_input_file();
        void copy_input_data();
        void alloc_base_attributes(int rank);
        void free_base_attributes();
};
#endif // NBODYSYSTEM_HPP
