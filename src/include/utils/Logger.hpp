#ifndef LOGGER_HPP
#define LOGGER_HPP
#include "../common.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

class Logger {
    public:
        Logger(int print_log, std::ofstream &out_file);
        ~Logger();

        std::ostream *gstream;

        void print_info(int n, double e2, double eta, float itime, float hmr_time,
                        float cr_time);
        void print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f,
                       double *dt);
        void print_energy_log(double ITIME, int iterations, long long interactions,
                              int nsteps, Gtime gtime, Energy &energy,
                              double new_energy);
        void print_lagrange_radii(double ITIME, std::vector<double> lagrange_radii);
};

#endif
