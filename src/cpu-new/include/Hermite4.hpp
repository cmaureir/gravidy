#ifndef HERMITE4_HPP
#define HERMITE4_HPP
#include "common.hpp"
#include <cmath>
#include <limits>
#include "extra_utils.hpp"

class Hermite4 {
    public:
        Hermite4(int n, double e2, float eta);
        ~Hermite4();

        int n;
        double e2;
        float eta;

        int  find_particles_to_move(int *move, double ITIME, double *dt, double *t);
        void next_integration_time(double &ATIME, double *dt, double *t);
        void init_dt(double &ATIME, double *dt, double *t, Forces *f);
        void save_old_acc_jrk(int nact, int *move, Forces *old, Forces *f);

        void force_calculation(int i, int j, Predictor *p, Forces *f);
        void init_acc_jrk(Predictor *p, Forces* f);
        void predicted_pos_vel(double ITIME, Predictor *p, double4 *r, double4 *v,
                              Forces *f, double *t, Gtime &gtime);
        void update_acc_jrk(int nact, int *move, Predictor *p, Forces* f, Gtime &gtime);
        void correction_pos_vel(double ITIME, int nact, int *move, double4 *r,
                                double4 *v, Forces *f, double *t, double *dt,
                                Predictor *p, Forces *old, double4 *a3, double4 *a2,
                                Gtime &gtime);

};

#endif
