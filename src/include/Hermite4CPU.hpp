#ifndef HERMITE4CPU_HPP
#define HERMITE4CPU_HPP
#include "Hermite4.hpp"

class Hermite4CPU : public Hermite4 {
    public:
        Hermite4CPU(int n, double e2, float eta) : Hermite4(n, e2, eta) { }

        void predicted_pos_vel(double ITIME, Predictor *p, double4 *r, double4 *v,
                               Forces *f, double *t, Gtime &gtime);
        void correction_pos_vel(double ITIME, int nact, int *move, double4 *r,
                                double4 *v, Forces *f, double *t, double *dt,
                                Predictor *p, Forces *old, double4 *a3, double4 *a2,
                                Gtime &gtime);
        void force_calculation(int i, int j, Predictor *p, Forces *f);
        void init_acc_jrk(Predictor *p, Forces* f);
        void update_acc_jrk(int nact, int *move, Predictor *p, Forces* f, Gtime &gtime);

};

#endif
