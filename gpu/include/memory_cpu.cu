#include "memory_cpu.hpp"

/*
 * @fn alloc_vectors_cpu()
 *
 */
void alloc_vectors_cpu()
{
    d4_size = n * sizeof(double4);
    d1_size = n * sizeof(double);
    f1_size = n * sizeof(float);
    i1_size = n * sizeof(int);
    size_t fsize = n * sizeof(Forces);
    size_t psize = n * sizeof(Predictor);

    /*
     * CPU pointers
     */
    cudaHostAlloc((void**)&h_r,        d4_size,         cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_v,        d4_size,         cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_f,        fsize,           cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_p,        psize,           cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_move,     i1_size,         cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_fout,     fsize * NJBLOCK, cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_fout_tmp, fsize * NJBLOCK, cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_i,        psize,           cudaHostAllocDefault);

    h_a2     = (double4*) malloc(d4_size);
    h_a3     = (double4*) malloc(d4_size);
    h_old_a  = (double4*) malloc(d4_size);
    h_old_a1 = (double4*) malloc(d4_size);
    h_ekin   =  (double*) malloc(d1_size);
    h_epot   =  (double*) malloc(d1_size);
    h_t      =  (double*) malloc(d1_size);
    h_dt     =  (double*) malloc(d1_size);
    h_m      =   (float*) malloc(f1_size);


    memset(h_f     , 0, fsize);
    memset(h_a2    , 0, d4_size);
    memset(h_a3    , 0, d4_size);
    memset(h_old_a , 0, d4_size);
    memset(h_old_a1, 0, d4_size);
    memset(h_t     , 0, d1_size);
    memset(h_dt    , 0, d1_size);
    memset(h_move  , 0, i1_size);

    for (int i = 0; i < n; i++)
    {
        h_m[i]      = part[i].m;
        h_r[i].x    = part[i].r.x;
        h_r[i].y    = part[i].r.y;
        h_r[i].z    = part[i].r.z;
        h_v[i].x    = part[i].v.x;
        h_v[i].y    = part[i].v.y;
        h_v[i].z    = part[i].v.z;
        h_p[i].r[0] = part[i].r.x;
        h_p[i].r[1] = part[i].r.y;
        h_p[i].r[2] = part[i].r.z;
        h_p[i].v[0] = part[i].v.x;
        h_p[i].v[1] = part[i].v.y;
        h_p[i].v[2] = part[i].v.z;
    }
}


/*
 * @fn clean_vectors_cpu()
 *
 * @brief
 *  Free memory on the CPU
 */
void free_vectors_cpu()
{
    cudaFreeHost(h_r);
    cudaFreeHost(h_v);
    cudaFreeHost(h_f);
    cudaFreeHost(h_p);
    cudaFreeHost(h_move);

    free(h_a2);
    free(h_a3);
    free(h_old_a);
    free(h_old_a1);
    free(h_ekin);
    free(h_epot);
    free(h_t);
    free(h_dt);
    free(h_m);

    cudaFreeHost(h_fout);
    cudaFreeHost(h_fout_tmp);
    cudaFreeHost(h_i);
}
