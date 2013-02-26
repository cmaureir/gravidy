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
    h_r      = (double4*) malloc(d4_size);
    h_v      = (double4*) malloc(d4_size);
    h_f      = (Forces*) malloc(sizeof(Forces) * n);
    h_a2     = (double4*) malloc(d4_size);
    h_a3     = (double4*) malloc(d4_size);
    h_old_a  = (double4*) malloc(d4_size);
    h_old_a1 = (double4*) malloc(d4_size);
    h_p      = (Predictor*) malloc(sizeof(Predictor) * n);
    h_ekin   =  (double*) malloc(d1_size);
    h_epot   =  (double*) malloc(d1_size);
    h_t      =  (double*) malloc(d1_size);
    h_dt     =  (double*) malloc(d1_size);
    h_m      =   (float*) malloc(f1_size);
    h_move   =     (int*) malloc(i1_size);

    // Empty double4
    double4 empty = {0.0, 0.0, 0.0, 0.0};

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

        h_f[i].a[0]  = 0.0;
        h_f[i].a[1]  = 0.0;
        h_f[i].a[2]  = 0.0;
        h_f[i].a1[0] = 0.0;
        h_f[i].a1[1] = 0.0;
        h_f[i].a1[2] = 0.0;
        h_a2[i]      = empty;
        h_a3[i]      = empty;
        h_old_a[i]   = empty;
        h_old_a1[i]  = empty;
        h_t[i]       = 0.0;
        h_dt[i]      = 0.0;
        h_move[i]    = 0;
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
    free(h_r);
    free(h_v);
    free(h_f);
    free(h_a2);
    free(h_a3);
    free(h_t);
    free(h_dt);
    free(h_old_a);
    free(h_old_a1);
    free(h_p);
    free(h_ekin);
    free(h_epot);
    free(h_move);
}
