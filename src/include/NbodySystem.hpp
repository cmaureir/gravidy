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
#include "utils/Logger.hpp"
#include "utils/NbodyUtils.hpp"

#ifdef GPU
#include "gpu/Hermite4GPU.cuh"
#elif KEPLER
#include "kepler/Hermite4Kepler.hpp"
#else
#include "Hermite4CPU.hpp"
#endif

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
        Predictor *h_i;
        double    *h_ekin;
        double    *h_epot;
        Forces    *h_fout_tmp;

        // Device Particles arrays
        //float     *d_m;
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


        void read_input_file();

        void alloc_arrays_host();
        void alloc_arrays_device();
        void free_arrays_host();
        void free_arrays_device();

        void copy_initial_data();
        //void integration(Hermite4 *h4, Logger log, NbodyUtils nu);

        /** Normal Hermite */
        // GPU
        #ifdef GPU
        void integration(Hermite4GPU h4, Logger log, NbodyUtils nu);
        double get_energy_gpu();
        /** Hermite/Kepler */
        #elif KEPLER
        void integration(Hermite4Kepler h4, Logger log, NbodyUtils nu);
        #else
        // CPU
        void integration(Hermite4CPU h4, Logger log, NbodyUtils nu);
        #endif
        double get_energy();

};

#endif
