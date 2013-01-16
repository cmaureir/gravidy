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

    /*
     * CPU pointers
     */
    //h_r      = (double4*) malloc(d4_size);
    //h_v      = (double4*) malloc(d4_size);
    //h_a      = (double4*) malloc(d4_size);
    //h_a1     = (double4*) malloc(d4_size);
    cudaHostAlloc((void**)&h_r,  d4_size, cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_v,  d4_size, cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_a,  d4_size, cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_a1, d4_size, cudaHostAllocDefault);
    h_a2     = (double4*) malloc(d4_size);
    h_a3     = (double4*) malloc(d4_size);
    h_old_a  = (double4*) malloc(d4_size);
    h_old_a1 = (double4*) malloc(d4_size);
    //h_p_r    = (double4*) malloc(d4_size);
    //h_p_v    = (double4*) malloc(d4_size);
    cudaHostAlloc((void**)&h_p_r,  d4_size, cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_p_v,  d4_size, cudaHostAllocDefault);
    h_ekin   =  (double*) malloc(d1_size);
    h_epot   =  (double*) malloc(d1_size);
    h_t      =  (double*) malloc(d1_size);
    h_dt     =  (double*) malloc(d1_size);
    h_m      =   (float*) malloc(f1_size);
    //h_move   =     (int*) malloc(i1_size);
    cudaHostAlloc((void**)&h_move,  i1_size, cudaHostAllocDefault);

    // Empty double4
    double4 empty = {0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < n; i++)
    {

        h_m[i]     = part[i].m;
        h_r[i].x   = part[i].r.x;
        h_r[i].y   = part[i].r.y;
        h_r[i].z   = part[i].r.z;
        h_v[i].x   = part[i].v.x;
        h_v[i].y   = part[i].v.y;
        h_v[i].z   = part[i].v.z;
        h_p_r[i].x = part[i].r.x;
        h_p_r[i].y = part[i].r.y;
        h_p_r[i].z = part[i].r.z;
        h_p_v[i].x = part[i].v.x;
        h_p_v[i].y = part[i].v.y;
        h_p_v[i].z = part[i].v.z;

        h_a[i]      = empty;
        h_a1[i]     = empty;
        h_a2[i]     = empty;
        h_a3[i]     = empty;
        h_old_a[i]  = empty;
        h_old_a1[i] = empty;
        h_t[i]      = 0.0;
        h_dt[i]     = 0.0;
        h_move[i]   = 0;
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
    free(h_m);
    //free(h_r);
    //free(h_v);
    //free(h_a);
    //free(h_a1);
    cudaFreeHost(h_r);
    cudaFreeHost(h_v);
    cudaFreeHost(h_a);
    cudaFreeHost(h_a1);
    free(h_a2);
    free(h_a3);
    free(h_t);
    free(h_dt);
    free(h_old_a);
    free(h_old_a1);
    //free(h_p_r);
    //free(h_p_v);
    cudaFreeHost(h_p_r);
    cudaFreeHost(h_p_v);
    free(h_ekin);
    free(h_epot);
    //free(h_move);
    cudaFreeHost(h_move);
}
