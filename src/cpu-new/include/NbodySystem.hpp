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
#include "OptionsParser.hpp"
#include "Hermite4.hpp"
#include "Logger.hpp"
#include "NbodyUtils.hpp"

class NbodySystem {
    public:
        NbodySystem(OptionsParser op);
        ~NbodySystem();

        Gtime gtime;
        std::vector<file_data> reader;
        int print_log;

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

        // Host Particles attributes
        //float     *h_m;
        double4   *h_r;
        double4   *h_v;
        Forces    *h_f;
        double4   *h_a2;
        double4   *h_a3;
        Forces    *h_old;
        double    *h_t;
        double    *h_dt;
        int       *h_move;
        Predictor *h_p;
        double    *h_ekin;
        double    *h_epot;

        void read_input_file();
        void alloc_arrays_host();
        void copy_initial_data();
        void integration(Hermite4 h4, Logger log, NbodyUtils nu);
        void free_arrays_host();
        double get_energy();
};

#endif
