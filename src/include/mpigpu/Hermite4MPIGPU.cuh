#ifndef HERMITE4MPI_HPP
#define HERMITE4MPI_HPP
#include "../Hermite4.hpp"

#define MPI_NUM_SLAVES 32

#define CUDA_SAFE_CALL_NO_SYNC( call) do {                          \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
    } } while (0)
#define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);

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

class Hermite4MPIGPU : public Hermite4 {
    public:
        Hermite4MPIGPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu,
                    int rank, int nproc);
        ~Hermite4MPIGPU();

        size_t nthreads;
        size_t nblocks;
        size_t smem;
        size_t smem_reduce;

        cudaEvent_t start;
        cudaEvent_t stop;

        int num_devices, num_cores;
        int rank;
        int nprocs;
        int tag;
        MPI_Status   status;
        MPI_Datatype f_type;
        MPI_Op       f_op;

        int chunk_size;
        int chunk_begin;
        int chunk_end;
        Forces *h_tmp_f;
        Forces *d_tmp_f;

        void alloc_slaves_memory(int rank);
        void force_calculation(Predictor pi, Predictor pj, Forces &fi);
        void init_acc_jrk();
        void update_acc_jrk(int nact);
        void predicted_pos_vel(double ITIME);
        void correction_pos_vel(double ITIME, int nact);
        void integration();
};

__global__ void k_init_acc_jrk(Predictor *p,
                               Forces *f,
                               int n,
                               double e2,
                               int chunk_begin,
                               int chunk_end);

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
