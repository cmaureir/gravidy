#include "memory_gpu.cuh"

/*
 * @fn alloc_vectors_gpu()
 *
 */
void alloc_vectors_gpu()
{
    /*
     * GPU pointers
     */
    size_t fsize = sizeof(Forces) * n;
    size_t psize = sizeof(Predictor) * n;

    CUDA_SAFE_CALL(cudaMalloc((void**)&d_r,     d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_v,     d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_f,     fsize));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_p,     psize));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_ekin,  d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_epot,  d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_t,     d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_dt,    d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_m,     f1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_move,  i1_size));

    CUDA_SAFE_CALL(cudaMalloc((void**)&d_i,    psize));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_fout, fsize * NJBLOCK));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_fout_tmp, fsize * NJBLOCK));
    CUDA_SAFE_CALL(cudaMemset(d_i,  0, psize));
    CUDA_SAFE_CALL(cudaMemset(d_fout,  0, fsize * NJBLOCK));
    CUDA_SAFE_CALL(cudaMemset(d_fout_tmp,  0, fsize * NJBLOCK));

    /*
     * Memset
     */
    CUDA_SAFE_CALL(cudaMemset(d_r, 0, d4_size));
    CUDA_SAFE_CALL(cudaMemset(d_v, 0, d4_size));
    CUDA_SAFE_CALL(cudaMemset(d_f, 0, fsize));
    CUDA_SAFE_CALL(cudaMemset(d_p, 0, psize));
    CUDA_SAFE_CALL(cudaMemset(d_ekin,  0,d1_size));
    CUDA_SAFE_CALL(cudaMemset(d_epot,  0,d1_size));
    CUDA_SAFE_CALL(cudaMemset(d_t,     0,d1_size));
    CUDA_SAFE_CALL(cudaMemset(d_dt,    0,d1_size));
    CUDA_SAFE_CALL(cudaMemset(d_m,     0,f1_size));
    CUDA_SAFE_CALL(cudaMemset(d_move,  0,i1_size));
}


/*
 * @fn free_vectors_gpu()
 *
 * @brief
 *  Free memory on the GPU
 */
void free_vectors_gpu()
{
    CUDA_SAFE_CALL(cudaFree(d_r));
    CUDA_SAFE_CALL(cudaFree(d_v));
    CUDA_SAFE_CALL(cudaFree(d_m));
    CUDA_SAFE_CALL(cudaFree(d_t));
    CUDA_SAFE_CALL(cudaFree(d_f));
    CUDA_SAFE_CALL(cudaFree(d_p));
    CUDA_SAFE_CALL(cudaFree(d_dt));
    CUDA_SAFE_CALL(cudaFree(d_ekin));
    CUDA_SAFE_CALL(cudaFree(d_epot));
    CUDA_SAFE_CALL(cudaFree(d_move));

    CUDA_SAFE_CALL(cudaFree(d_i));
    CUDA_SAFE_CALL(cudaFree(d_fout));
}
