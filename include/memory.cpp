#include "memory.hpp"

/*
 * @fn init_vectors()
 *
 * @brief
 *  Memory allocation in host (CPU) and device (GPU)
 */
void init_vectors()
{
    int d4_size = n * sizeof(double4);
    int f1_size = n * sizeof(float);
    int i1_size = n * sizeof(int);

    h_r     = (double4*) malloc(d4_size);
    h_v     = (double4*) malloc(d4_size);
    h_a     = (double4*) malloc(d4_size);
    h_j     = (double4*) malloc(d4_size);
    h_old_a = (double4*) malloc(d4_size);
    h_old_j = (double4*) malloc(d4_size);
    h_new_a = (double4*) malloc(d4_size);
    h_new_j = (double4*) malloc(d4_size);
    h_p_r   = (double4*) malloc(d4_size);
    h_p_v   = (double4*) malloc(d4_size);
    h_ekin  =  (float*) malloc(f1_size);
    h_epot  =  (float*) malloc(f1_size);
    h_t     =  (float*) malloc(f1_size);
    h_dt    =  (float*) malloc(f1_size);
    h_m     =   (float*) malloc(f1_size);
    h_move  =     (int*) malloc(i1_size);

    cudaMalloc((void**)&d_r,     d4_size);
    cudaMalloc((void**)&d_v,     d4_size);
    cudaMalloc((void**)&d_a,     d4_size);
    cudaMalloc((void**)&d_j,     d4_size);
    cudaMalloc((void**)&d_old_a, d4_size);
    cudaMalloc((void**)&d_old_j, d4_size);
    cudaMalloc((void**)&d_new_a, d4_size);
    cudaMalloc((void**)&d_new_j, d4_size);
    cudaMalloc((void**)&d_p_r,   d4_size);
    cudaMalloc((void**)&d_p_v,   d4_size);
    cudaMalloc((void**)&d_m,     f1_size);
    cudaMalloc((void**)&d_ekin,  f1_size);
    cudaMalloc((void**)&d_epot,  f1_size);
    cudaMalloc((void**)&d_t,     f1_size);
    cudaMalloc((void**)&d_dt,    f1_size);
    cudaMalloc((void**)&d_move,  i1_size);

    for (int i = 0; i < n; i++)
    {
        h_r[i].x = part[i].r.x;
        h_r[i].y = part[i].r.y;
        h_r[i].z = part[i].r.z;

        h_v[i].x = part[i].v.x;
        h_v[i].y = part[i].v.y;
        h_v[i].z = part[i].v.z;

        h_a[i].x = 0;
        h_a[i].y = 0;
        h_a[i].z = 0;

        h_j[i].x = 0;
        h_j[i].y = 0;
        h_j[i].z = 0;

        h_old_a[i].x = 0;
        h_old_a[i].y = 0;
        h_old_a[i].z = 0;

        h_new_a[i].x = 0;
        h_new_a[i].y = 0;
        h_new_a[i].z = 0;

        h_old_j[i].x = 0;
        h_old_j[i].y = 0;
        h_old_j[i].z = 0;

        h_new_j[i].x = 0;
        h_new_j[i].y = 0;
        h_new_j[i].z = 0;

        h_p_r[i].x = 0;
        h_p_r[i].y = 0;
        h_p_r[i].z = 0;

        h_p_v[i].x = 0;
        h_p_v[i].y = 0;
        h_p_v[i].z = 0;

        h_t[i]  = 0.0f;
        h_dt[i] = 0.0f;

        h_move[i] = 0;

        h_m[i]   = part[i].m;

    }

    // The mass is always the same, so to avoid copying it every
    //  function, we copy it at the begining.

    cudaMemcpy(d_m, h_m, f1_size, cudaMemcpyHostToDevice);
}


/*
 * @fn clean_vectors()
 *
 * @brief
 *  Free memory on the GPU and CPU
 */
void clean_vectors()
{

    free(h_r);
    free(h_v);
    free(h_a);
    free(h_j);
    free(h_m);
    free(h_t);
    free(h_dt);
    free(h_old_a);
    free(h_old_j);
    free(h_new_a);
    free(h_new_j);
    free(h_p_r);
    free(h_p_v);
    free(h_ekin);
    free(h_epot);
    free(h_move);

    cudaFree(d_r);
    cudaFree(d_v);
    cudaFree(d_a);
    cudaFree(d_j);
    cudaFree(d_m);
    cudaFree(d_t);
    cudaFree(d_old_a);
    cudaFree(d_old_j);
    cudaFree(d_new_a);
    cudaFree(d_new_j);
    cudaFree(d_p_r);
    cudaFree(d_p_v);
    cudaFree(d_dt);
    cudaFree(d_ekin);
    cudaFree(d_epot);
    cudaFree(d_move);
}
