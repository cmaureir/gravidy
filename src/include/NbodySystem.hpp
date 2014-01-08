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
        int n;
        long long int iterations;
        double e2;
        float eta;
        float total_mass;
        float integration_time;

        // Global parameters
        Energy en;
        float hmr_time;
        float cr_time;

        // Files
        std::string input_filename;
        std::string output_filename;
        std::ofstream out_file;

        // Host Particles arrays
        double4   *h_r;
        double4   *h_v;
        Forces    *h_f;
        double4   *h_a2;
        double4   *h_a3;
        Forces    *h_old;
        double    *h_t;
        double    *h_dt;
        double    *h_dt_old;
        int       *h_move;
        Predictor *h_p;
        Predictor *h_i;
        double    *h_ekin;
        double    *h_epot;
        Forces    *h_fout_tmp;

        // Device Particles arrays
        double4   *d_r;
        double4   *d_v;
        Forces    *d_f;
        Forces    *d_fout;
        Forces    *d_fout_tmp;
        Forces    *d_old;
        double    *d_t;
        double    *d_dt;
        int       *d_move;
        Predictor *d_p;
        Predictor *d_i;
        double    *d_ekin;
        double    *d_epot;

        // For keplerian correction
        Forces    *h_fbh;
        double4   *h_a2bh;
        double4   *h_a3bh;

        void read_input_file();
        void alloc_base_attributes();
        void free_base_attributes();
        void copy_input_data();
};
#endif // NBODYSYSTEM_HPP
