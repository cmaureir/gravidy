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

class Hermite4GPU : public Hermite4 {
    public:
        Hermite4GPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
            : Hermite4(ns, logger, nu)
        {
            smem = sizeof(Predictor) * BSIZE;
            smem_reduce = sizeof(Forces) * NJBLOCK + 1;
            gpus = 0;

            CSC(cudaGetDeviceCount(&gpus));
            std::cout << "GPUs: " << gpus << std::endl;
            gpus = 1;
            // cudaSetDevice(ID); // set the active gpu
            // cudaGetDevice(&ID) // get the current gpu

            std::cout << "Spliting " << ns->n << " particles in " << gpus << " GPUs" << std::endl;

            if (ns->n % gpus == 0)
            {
                int size = ns->n/gpus;
                for ( int g = 0; g < gpus; g++)
                    n_part[g] = size;
            }
            else
            {
                int size = std::ceil(ns->n/(float)gpus);
                for ( int g = 0; g < gpus; g++)
                {
                    if (ns->n - size*(g+1) > 0)
                        n_part[g] = size;
                    else
                        n_part[g] = ns->n - size*g;
                }
            }

            for(int g = 0; g < gpus; g++)
            {
                std::cout << "GPU " << g << " particles: " << n_part[g] << std::endl;
            }

            alloc_arrays_device();
        }
        ~Hermite4GPU();
        using Hermite4::init_acc_jrk;
        using Hermite4::update_acc_jrk;
        using Hermite4::predicted_pos_vel;
        using Hermite4::correction_pos_vel;

        size_t nthreads;
        size_t nblocks;
        size_t smem;
        size_t smem_reduce;

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

__global__ void k_init_acc_jrk(Predictor *p,
                               Forces *f,
                               int n,
                               double e2,
                               int dev,
                               int dev_size);

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
