#ifndef HERMITE4_HPP
#define HERMITE4_HPP
#include "common.hpp"
#include <cmath>
#include "extra_utils.hpp"

class Hermite4 {
    public:
        Hermite4(int n, float e2);
        ~Hermite4();

        int n;
        float e2;

        void init_acc_jrk(double4* r, double4* v, Forces* f);
        void force_calculation(int i, int j, double4 *r, double4 *v, Forces *f);
        void init_dt(double &ATIME, Forces *f, double *dt, double *t);
        int find_particles_to_move(double ITIME, double *dt, double *t, int *move);
        void save_old_acc_jrk(Forces *f, int *move, int nact, float4 *old_a, float4 *old_a1);
        void predicted_pos_vel(double ITIME, Predictor *p, double4 *r, double4 *v, Forces *f, double *t);
        void update_acc_jrk(int nact, int *move, double4* r, double4* v, Forces* f);
        void correction_pos_vel(double ITIME, int nact, int *move, double4 *r, double4 *v, Forces *f, double *t, double *dt, Predictor *p, float4 *old_a, float4 *old_a1, float4 *a3, float4 *a2, float eta);

};

#endif
