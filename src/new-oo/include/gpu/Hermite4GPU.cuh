#ifndef HERMITE4GPU_HPP
#define HERMITE4GPU_HPP
#include "../Hermite4.hpp"
//#include <string>

class Hermite4GPU : public Hermite4 {
    public:
        Hermite4GPU(int n, double e2, float eta) : Hermite4(n, e2, eta) {
            nthreads = BSIZE;
            nblocks = std::ceil(n/(float)nthreads);
            }

        size_t nthreads, nblocks;
        cudaEvent_t start, stop;

        Predictor *d_p;
        Predictor *d_i;
        Predictor *h_i;
        Forces *d_fout;
        Forces *d_fout_tmp;
        Forces *h_fout_tmp;
        Forces *d_f;
        int *d_move;

        void get_kernel_error();
        void gpu_timer_start();
        float gpu_timer_stop(std::string f);

        void set_pointers(Predictor*, Predictor*, Predictor*, Forces*, Forces*,
                          Forces*, Forces*, int*);
        void init_acc_jrk(Predictor *p, Forces* f);
        void update_acc_jrk(int nact, int *move, Predictor *p, Forces* f, Gtime &gtime);
};

__global__ void k_energy(double4*,
                         double4*,
                         double*,
                         double*,
                         int,
                         double);

__global__ void k_init_acc_jrk(Predictor*,
                               Forces*,
                               int,
                               double);

__device__ void k_force_calculation(Predictor,
                                     Predictor,
                                     Forces&,
                                     double);

__global__ void k_update(Predictor*,
                         Predictor*,
                         Forces*,
                         int*,
                         int,
                         int,
                         double);

__global__ void reduce(Forces*,
                       Forces*,
                       unsigned int);

#endif
