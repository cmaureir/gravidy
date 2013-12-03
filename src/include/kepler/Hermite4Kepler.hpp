#ifndef HERMITE4KEPLER_HPP
#define HERMITE4KEPLER_HPP
#include "../Hermite4.hpp"

#include <cmath>
#include <iostream>

static double
    _1_17_16 = 1./17./16.,
    _1_16_15 = 1./16./15.,
    _1_15_14 = 1./15./14.,
    _1_14_13 = 1./14./13.,
    _1_13_12 = 1./13./12.,
    _1_12_11 = 1./12./11.,
    _1_11_10 = 1./11./10.,
    _1_10_9  = 1./10./9.,
    _1_9_8   = 1./9./8.,
    _1_8_7   = 1./8./7.,
    _1_7_6   = 1./7./6.,
    _1_6_5   = 1./6./5.,
    _1_4_3   = 1./4./3.,
    _1_3_2   = 1./3./2.;

class Hermite4Kepler : public Hermite4 {
    public:
        Hermite4Kepler(int n, double e2, float eta) : Hermite4(n, e2, eta) { }

        // Hermite
        void predicted_pos_vel(double ITIME, Predictor *p, double4 *r, double4 *v,
                               Forces *f, double *t, Gtime &gtime);
        void correction_pos_vel(double ITIME, int nact, int *move, double4 *r,
                                double4 *v, Forces *f, double *t, double *dt,
                                Predictor *p, Forces *old, double4 *a3, double4 *a2,
                                Gtime &gtime);
        void force_calculation(int i, int j, Predictor *p, Forces *f);
        void init_acc_jrk(Predictor *p, Forces* f);
        void update_acc_jrk(int nact, int *move, Predictor *p, Forces* f, Gtime &gtime);

        // Kepler
        void predicted_pos_vel_kepler(double, int);
        void kepler_prediction(double*, double*, double*, double*,
                               double*, double*, double,  int);
        double solve_kepler(double, double);
        double kepler(const double, double);
};

#endif
