#ifndef NBODYSYSTEM_HPP
#define NBODYSYSTEM_HPP
#include <vector>
#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <algorithm>

#include "common.hpp"
#include "OptionsParser.hpp"
#include "Hermite4.hpp"

class NbodySystem {
    public:
        NbodySystem();
        ~NbodySystem();

        // Misc
        Gtime gtime;
        std::vector<file_data> reader;
        int print_log;

        // Configuration parameters
        int n;
        int iterations;
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
        FILE *out;

        // Host Particles attributes
        //float     *h_m;
        double4   *h_r;
        double4   *h_v;
        Forces    *h_f;
        double4   *h_a2;
        double4   *h_a3;
        double4   *h_old_a;
        double4   *h_old_a1;
        double    *h_t;
        double    *h_dt;
        int       *h_move;
        Predictor *h_p;
        double    *h_ekin;
        double    *h_epot;

        void get_parameters(OptionsParser op);
        void read_input_file();
        void alloc_arrays_host();
        void copy_initial_data();
        void free_arrays_host();
        void integration(Hermite4 h4);
        double get_energy();
};

#endif
