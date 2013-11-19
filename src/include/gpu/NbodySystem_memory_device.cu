#include "../NbodySystem.hpp"

void NbodySystem::alloc_arrays_device()
{
    int d4_size = n * sizeof(double4);
    int d1_size = n * sizeof(double);
    int i1_size = n * sizeof(int);
    int ff_size = n * sizeof(Forces);
    int pp_size = n * sizeof(Predictor);

    CUDA_SAFE_CALL(cudaMalloc((void**)&d_r,        d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_v,        d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_f,        ff_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_p,        pp_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_ekin,     d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_epot,     d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_t,        d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_dt,       d1_size));
    //CUDA_SAFE_CALL(cudaMalloc((void**)&d_m,        f1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_move,     i1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_i,        pp_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_fout,     ff_size * NJBLOCK));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_fout_tmp, ff_size * NJBLOCK));

    CUDA_SAFE_CALL(cudaMemset(d_r,         0, d4_size));
    CUDA_SAFE_CALL(cudaMemset(d_v,         0, d4_size));
    CUDA_SAFE_CALL(cudaMemset(d_f,         0, ff_size));
    CUDA_SAFE_CALL(cudaMemset(d_p,         0, pp_size));
    CUDA_SAFE_CALL(cudaMemset(d_ekin,      0, d1_size));
    CUDA_SAFE_CALL(cudaMemset(d_epot,      0, d1_size));
    CUDA_SAFE_CALL(cudaMemset(d_t,         0, d1_size));
    CUDA_SAFE_CALL(cudaMemset(d_dt,        0, d1_size));
    //CUDA_SAFE_CALL(cudaMemset(d_m,         0, f1_size));
    CUDA_SAFE_CALL(cudaMemset(d_move,      0, i1_size));
    CUDA_SAFE_CALL(cudaMemset(d_i,         0, pp_size));
    CUDA_SAFE_CALL(cudaMemset(d_fout,      0, ff_size * NJBLOCK));
    CUDA_SAFE_CALL(cudaMemset(d_fout_tmp,  0, ff_size * NJBLOCK));

}

void NbodySystem::free_arrays_device()
{
    CUDA_SAFE_CALL(cudaFree(d_r));
    CUDA_SAFE_CALL(cudaFree(d_v));
    CUDA_SAFE_CALL(cudaFree(d_f));
    CUDA_SAFE_CALL(cudaFree(d_p));
    CUDA_SAFE_CALL(cudaFree(d_ekin));
    CUDA_SAFE_CALL(cudaFree(d_epot));
    CUDA_SAFE_CALL(cudaFree(d_t));
    CUDA_SAFE_CALL(cudaFree(d_dt));
    //CUDA_SAFE_CALL(cudaFree(d_m));
    CUDA_SAFE_CALL(cudaFree(d_move));
    CUDA_SAFE_CALL(cudaFree(d_i));
    CUDA_SAFE_CALL(cudaFree(d_fout));
    CUDA_SAFE_CALL(cudaFree(d_fout_tmp));
}
