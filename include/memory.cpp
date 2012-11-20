#include "memory.hpp"

/*
 * @fn init_vectors()
 *
 */
void init_vectors()
{
    d4_size = n * sizeof(double4);
    d1_size = n * sizeof(double);
    f1_size = n * sizeof(float);
    i1_size = n * sizeof(int);

    h_r     = (double4*) malloc(d4_size);
    h_v     = (double4*) malloc(d4_size);

    h_a     = (double4*) malloc(d4_size);
    h_a1     = (double4*) malloc(d4_size);

    h_a2    = (double4*) malloc(d4_size);
    h_a3    = (double4*) malloc(d4_size);

    h_old_a = (double4*) malloc(d4_size);
    h_old_a1 = (double4*) malloc(d4_size);


    //h_new_a = (double4*) malloc(n * NJBLOCK * sizeof(double4));
    //h_new_j = (double4*) malloc(n * NJBLOCK * sizeof(double4));
    //h_new_a = (double4*) malloc(d4_size);
    //h_new_j = (double4*) malloc(d4_size);

    //h_new_a = (double4**) malloc(n * sizeof(double4*));
    //h_new_j = (double4**) malloc(n * sizeof(double4*));
    //for (int ii = 0; ii < n; ii++) {
    //    h_new_a[ii] = (double4*)malloc(NJBLOCK*sizeof(double4));
    //    h_new_j[ii] = (double4*)malloc(NJBLOCK*sizeof(double4));
    //}

    h_p_r   = (double4*) malloc(d4_size);
    h_p_v   = (double4*) malloc(d4_size);

    h_ekin  =  (double*) malloc(d1_size);
    h_epot  =  (double*) malloc(d1_size);

    h_t     =  (double*) malloc(d1_size);
    h_dt    =  (double*) malloc(d1_size);

    h_m     =   (float*) malloc(f1_size);
    h_move  =     (int*) malloc(i1_size);

    //d_new_a = new double4* [n];
    //d_new_j = new double4* [n];
    //for (int ii = 0; ii < n; ii++) {
    //    cudaMalloc((void**)&d_new_a[ii], NJBLOCK * sizeof(double4));
    //    cudaMalloc((void**)&d_new_j[ii], NJBLOCK * sizeof(double4));
    //}
    //cudaMalloc((void**)&d_new_a, NJBLOCK * n * sizeof(double4));
    //cudaMalloc((void**)&d_new_j, NJBLOCK * n * sizeof(double4));

    //CUDA_SAFE_CALL(cudaMalloc((void**)&d_new_a, d4_size));
    //CUDA_SAFE_CALL(cudaMalloc((void**)&d_new_j, d4_size));


    CUDA_SAFE_CALL(cudaMalloc((void**)&d_r,     d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_v,     d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_a,     d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_a1,     d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_old_a, d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_old_a1, d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_p_r,   d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_p_v,   d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_m,     f1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_ekin,  d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_epot,  d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_t,     d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_dt,    d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_move,  i1_size));
    //CUDA_SAFE_CALL(cudaMalloc((void**)&tmp_red, d4_size));

    double4 empty = {0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < n; i++)
    {
        h_r[i].x = part[i].r.x;
        h_r[i].y = part[i].r.y;
        h_r[i].z = part[i].r.z;

        h_v[i].x = part[i].v.x;
        h_v[i].y = part[i].v.y;
        h_v[i].z = part[i].v.z;

        // Initial predicted position and velocity
        // has the same value, to perform
        // the energy calculation
        h_p_r[i].x = part[i].r.x;
        h_p_r[i].y = part[i].r.y;
        h_p_r[i].z = part[i].r.z;

        h_p_v[i].x = part[i].v.x;
        h_p_v[i].y = part[i].v.y;
        h_p_v[i].z = part[i].v.z;

        h_a[i] = empty;
        h_a1[i] = empty;

        h_a2[i] = empty;
        h_a3[i] = empty;

        h_old_a[i] = empty;
        h_old_a1[i] = empty;
        //for (int ii = 0; ii < NJBLOCK * n; ii++) {
        //    h_new_a[ii] = empty;
        //    h_new_j[ii] = empty;
        //}
        //h_new_a[i] = empty;
        //h_new_j[i] = empty;

        h_t[i]  = 0.0;
        h_dt[i] = 0.0;
        h_move[i] = 0;
        h_m[i]   = part[i].m;
    }
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
    free(h_a1);
    free(h_a2);
    free(h_a3);
    free(h_m);
    free(h_t);
    free(h_dt);
    free(h_old_a);
    free(h_old_a1);
    free(h_p_r);
    free(h_p_v);
    free(h_ekin);
    free(h_epot);
    free(h_move);

    CUDA_SAFE_CALL(cudaFree(d_r));
    CUDA_SAFE_CALL(cudaFree(d_v));
    CUDA_SAFE_CALL(cudaFree(d_a));
    CUDA_SAFE_CALL(cudaFree(d_a1));
    CUDA_SAFE_CALL(cudaFree(d_m));
    CUDA_SAFE_CALL(cudaFree(d_t));
    CUDA_SAFE_CALL(cudaFree(d_old_a));
    CUDA_SAFE_CALL(cudaFree(d_old_a1));
    CUDA_SAFE_CALL(cudaFree(d_p_r));
    CUDA_SAFE_CALL(cudaFree(d_p_v));
    CUDA_SAFE_CALL(cudaFree(d_dt));
    CUDA_SAFE_CALL(cudaFree(d_ekin));
    CUDA_SAFE_CALL(cudaFree(d_epot));
    CUDA_SAFE_CALL(cudaFree(d_move));

    //for (int i = 0; i < n; i++) {
    //    cudaFree(d_new_a[i]);
    //    cudaFree(d_new_j[i]);
    //    free(h_new_a[i]);
    //    free(h_new_j[i]);
    //}
    //CUDA_SAFE_CALL(cudaFree(d_new_a));
    //CUDA_SAFE_CALL(cudaFree(d_new_j));
    //CUDA_SAFE_CALL(cudaFree(tmp_red));
    //free(h_new_a);
    //free(h_new_j);
}
