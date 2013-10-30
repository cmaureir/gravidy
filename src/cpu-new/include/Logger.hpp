#ifndef LOGGER_HPP
#define LOGGER_HPP
#include "common.hpp"
#include <cstdio>

class Logger {
    public:
        Logger(int print_log);
        ~Logger();

        void print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f, double *dt);
};

#endif
