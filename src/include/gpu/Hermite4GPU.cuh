#ifndef HERMITE4GPU_HPP
#define HERMITE4GPU_HPP
#include "../Hermite4.hpp"
//#include <string>
#include <cassert>
#include <string>

class Hermite4GPU : public Hermite4 {
    public:
        Hermite4GPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
            : Hermite4(ns, logger, nu)
        {
            nthreads = BSIZE;
            nblocks = std::ceil(ns->n/(float)nthreads);
            smem = sizeof(Predictor) * BSIZE;
            smem_reduce = sizeof(Forces) * NJBLOCK + 1;

            alloc_arrays_device();
        }
        ~Hermite4GPU();

        size_t nthreads;
        size_t nblocks;
        size_t smem;
        size_t smem_reduce;

        cudaEvent_t start;
        cudaEvent_t stop;

        void alloc_arrays_device();
        void free_arrays_device();

        void force_calculation(int i, int j);
        void init_acc_jrk();
        void update_acc_jrk(int nact);
        void predicted_pos_vel(double ITIME);
        void correction_pos_vel(double ITIME, int nact);
        void integration();

        void get_kernel_error();
        void gpu_timer_start();
        float  gpu_timer_stop(std::string f);

        double get_energy_gpu();

};

__global__ void k_init_acc_jrk(Predictor *p,
                               Forces *f,
                               int n,
                               double e2);

__device__ void k_force_calculation(Predictor i_p,
                                    Predictor j_p,
                                    Forces &f,
                                    double e2);

__global__ void k_update(Predictor *i_p,
                         Predictor *j_p,
                         Forces *fout,
                         int *move,
                         int n,
                         int total,
                         double e2);

__global__ void reduce(Forces *in,
                       Forces *out);

__global__ void k_energy(double4 *r,
                         double4 *v,
                         double *ekin,
                         double *epot,
                         int n,
                         double e2);
#endif
