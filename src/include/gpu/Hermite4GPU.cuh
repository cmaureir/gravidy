/*
 * Copyright (c) 2016
 *
 * Cristián Maureira-Fredes <cmaureirafredes@gmail.com>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef HERMITE4GPU_HPP
#define HERMITE4GPU_HPP
#include "../Hermite4.hpp"
//#include <string>
#include <cassert>
#include <string>

#define CSC_NO_SYNC( call) do {                          \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
    } } while (0)
#define CSC( call)     CSC_NO_SYNC(call);

inline __host__ __device__ double4 operator+(const double4 &a, const double4 &b)
{
    double4 tmp = {a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w};
    return tmp;
}

inline __host__ __device__ void operator+=(double4 &a, double4 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
}

inline __host__ __device__ Forces operator+(Forces &a, Forces &b)
{
    Forces tmp;
    tmp.a[0] = a.a[0] + b.a[0];
    tmp.a[1] = a.a[1] + b.a[1];
    tmp.a[2] = a.a[2] + b.a[2];

    tmp.a1[0] = a.a1[0] + b.a1[0];
    tmp.a1[1] = a.a1[1] + b.a1[1];
    tmp.a1[2] = a.a1[2] + b.a1[2];

    return tmp;
}

inline __host__ __device__ void operator+=(Forces &a, Forces &b)
{
    a.a[0] += b.a[0];
    a.a[1] += b.a[1];
    a.a[2] += b.a[2];

    a.a1[0] += b.a1[0];
    a.a1[1] += b.a1[1];
    a.a1[2] += b.a1[2];
}

/**
 * Class which implements on the GPU the structure of the Hermite4 scheme.
 *
 * This contains all the implementations of the requirements processes to perform
 * the integration, like the initialization of the forces, the prediction,
 * the correction, and the general integration of the system.
 *
 */

class Hermite4GPU : public Hermite4 {
    public:
        Hermite4GPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu);
        ~Hermite4GPU();

        using Hermite4::init_acc_jrk;
        using Hermite4::update_acc_jrk;
        using Hermite4::predicted_pos_vel;
        using Hermite4::correction_pos_vel;
        using Hermite4::force_calculation;

        size_t nthreads;
        size_t nblocks;
        size_t smem;
        size_t smem_reduce;

        size_t i1_size;
        size_t d1_size;
        size_t d4_size;
        size_t ff_size;
        size_t pp_size;

        int gpus;
        int n_part[MAXGPUS];

        cudaEvent_t start;
        cudaEvent_t stop;

        void alloc_arrays_device();
        void free_arrays_device();

        void force_calculation(Predictor pi, Predictor pj, Forces &fi);
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

/**
 * Initialization kernel, which consider an \f$N^2\f$ interaction of the particles.
 */
__global__ void k_init_acc_jrk(Predictor *p,
                               Forces *f,
                               int n,
                               double e2,
                               int dev,
                               int dev_size);

/**
 * Predictor kernel, in charge of performing the prediction step of all the particles
 * on each integration step.
 */
__global__ void k_prediction(Forces *f,
                             double4 *r,
                             double4 *v,
                             double *t,
                             Predictor *p,
                             int dev_size,
                             double ITIME);

/**
 * Force interaction kernel, in charge of performing gravitational interaction
 * computation between two particles.
 */
__device__ void k_force_calculation(Predictor i_p,
                                    Predictor j_p,
                                    Forces &f,
                                    double e2);

/**
 * Force kernel, in charge of performing distribution of how the \f$N, N_{act}\f$
 * particles will be distributed among the GPUs.
 * This kernel calls the k_prediction kernel.
 */
__global__ void k_update(Predictor *i_p,
                         Predictor *j_p,
                         Forces *fout,
                         int n,
                         int total,
                         double e2);

/**
 * Force reduction kernel, in charge of summing up all the preliminary results
 * of the forces for the \f$N_{act}\f$ particles.
 */
__global__ void k_reduce(Forces *in,
                       Forces *out,
                       int shift_id,
                       int shift);

/**
 * Energy kernel, in charge of the calculation of the kinetic and potential
 * energy on the GPUs.
 */
__global__ void k_energy(double4 *r,
                         double4 *v,
                         double *ekin,
                         double *epot,
                         int n,
                         int dev_size,
                         int dev);
#endif
