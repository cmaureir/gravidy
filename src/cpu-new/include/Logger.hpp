#ifndef LOGGER_HPP
#define LOGGER_HPP
#include "common.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>

class Logger {
    public:
        Logger(int print_log, std::ofstream &out_file);
        ~Logger();

        std::ostream *gstream;

        void print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f,
                       double *dt);
        void print_energy_log(double ITIME, int iterations, long long interactions,
                              int nsteps, Gtime gtime, Energy &energy,
                              double new_energy);
};

#endif
