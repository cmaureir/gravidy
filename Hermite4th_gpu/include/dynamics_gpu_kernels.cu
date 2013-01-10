#include "dynamics_gpu_kernels.cuh"

__global__ void
k_energy(double4 *r, double4 *v, double *ekin, double *epot, float *m, int n)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j;
    double ekin_tmp = 0.0;
    double epot_tmp = 0.0;

    if (i < n)
    {
        epot_tmp = 0.0;
        for (j = i+1; j < n; j++)
        {
            double rx = r[j].x - r[i].x;
            double ry = r[j].y - r[i].y;
            double rz = r[j].z - r[i].z;
            double r2 = rx*rx + ry*ry + rz*rz;
            epot_tmp -= (m[i] * m[j]) * rsqrt(r2);
        }

        double vx = v[i].x * v[i].x;
        double vy = v[i].y * v[i].y;
        double vz = v[i].z * v[i].z;
        double v2 = vx + vy + vz;
        ekin_tmp = 0.5 * m[i] * v2;

        ekin[i] = ekin_tmp;
        epot[i] = epot_tmp;
    }
}

__global__ void k_update_acc_jrk_simple
(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int *move, int n, int total)
{

    double4 aa  = {0.0, 0.0, 0.0, 0.0};
    double4 jj  = {0.0, 0.0, 0.0, 0.0};
    double4 pos = {0.0, 0.0, 0.0, 0.0};
    double4 vel = {0.0, 0.0, 0.0, 0.0};
    float mj;

    int id = threadIdx.x + blockDim.x * blockIdx.x;

    if (id < total)
    {
        int id_move = move[id];
        int k;
        pos = r[id_move];
        vel = v[id_move];
        mj  = m[id_move];
        for (k = INIT_PARTICLE; k < n; k++)
        {
            if(id == k) continue;
            k_force_calculation(pos, vel, r[k], v[k], aa, jj, mj);
        }

     a[id_move] = aa;
     j[id_move] = jj;
    }
}
// ***
//__global__ void k_update_acc_jrk_simple
//(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int *move, int n, int total)
//{
//    extern __shared__ double4 sh[];
//    double4 *sr = (double4*)sh;
//    double4 *sv = (double4*)&sr[blockDim.x];
//
//    double4 aa  = {0.0, 0.0, 0.0, 0.0};
//    double4 jj  = {0.0, 0.0, 0.0, 0.0};
//    double4 pos = {0.0, 0.0, 0.0, 0.0};
//    double4 vel = {0.0, 0.0, 0.0, 0.0};
//    float mj;
//
//    int id = threadIdx.x + blockDim.x * blockIdx.x;
//    int tx = threadIdx.x;
//
//    int id_move = move[id];
//
//    //if(id_move != -1)
//    //{
//        pos = r[id_move];
//        vel = v[id_move];
//    //}
//
//    int tile = 0;
//    for (int i = 0; i < n; i += BSIZE)
//    {
//        int idx = tile * BSIZE + tx;
//
//        sr[tx]   = r[idx];
//        sv[tx]   = v[idx];
//        mj = m[idx];
//        __syncthreads();
//
//        for (int k = 0; k < BSIZE; k++)
//        {
//            if(id_move != -1)
//            {
//                k_force_calculation(pos, vel, sr[k], sv[k], aa, jj, mj);
//            }
//        }
//        __syncthreads();
//        tile++;
//    }
//
//    if(id_move != -1)
//    {
//        a[id_move] = aa;
//        j[id_move] = jj;
//    }
//}
//
/*
 * @fn k_init_acc_jrk
 *
 * @param to do
 *
 * @desc GPU Kernel to perform the calculation of the initial acceleration
 *       and jerk of the system.
 *
 * @note Working properly
 */
__global__ void k_init_acc_jrk
(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int n)
{

    extern __shared__ double4 sh[];
    double4 *sr = (double4*)sh;
    double4 *sv = (double4*)&sr[blockDim.x];

    double4 aa = {0.0, 0.0, 0.0, 0.0};
    double4 jj = {0.0, 0.0, 0.0, 0.0};
    float mj;

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int tx = threadIdx.x;

    if (id < n)
    {
        double4 pos = r[id];
        double4 vel = v[id];

        int tile = 0;
        for (int i = 0; i < n; i += BSIZE)
        {
            int idx = tile * BSIZE + tx;

            sr[tx]   = r[idx];
            sv[tx]   = v[idx];
            mj = m[idx];
            __syncthreads();

            for (int k = 0; k < BSIZE; k++)
            {
                k_force_calculation(pos, vel, sr[k], sv[k], aa, jj, mj);
            }
            __syncthreads();
            tile++;
        }

        a[id] = aa;
        j[id] = jj;
    }
}

/*
 * @fn k_force_calculation
 *
 * @desc GPU Kernel which calculates the interaction between
 *       a i-particle and a j-particle.
 *
 * @note Working properly.
 *
 */
__device__ void k_force_calculation(double4 i_pos, double4 i_vel,
                                    double4 j_pos, double4 j_vel,
                                    double4 &acc,  double4 &jrk,
                                    float   j_mass)
{
    double rx = j_pos.x - i_pos.x;
    double ry = j_pos.y - i_pos.y;
    double rz = j_pos.z - i_pos.z;

    double vx = j_vel.x - i_vel.x;
    double vy = j_vel.y - i_vel.y;
    double vz = j_vel.z - i_vel.z;

    double r2 = rx*rx + ry*ry + rz*rz + E2;
    double rinv = rsqrt(r2);
    double r2inv = rinv  * rinv;
    double r3inv = r2inv * rinv;
    double r5inv = r2inv * r3inv;
    double mr3inv = r3inv * j_mass;
    double mr5inv = r5inv * j_mass;

    double rv = rx*vx + ry*vy + rz*vz;

    acc.x += (rx * mr3inv);
    acc.y += (ry * mr3inv);
    acc.z += (rz * mr3inv);

    jrk.x += (vx * mr3inv - (3 * rv) * rx * mr5inv);
    jrk.y += (vy * mr3inv - (3 * rv) * ry * mr5inv);
    jrk.z += (vz * mr3inv - (3 * rv) * rz * mr5inv);
}

/*
 * @fn k_predicted_pos_vel
 *
 * @param to do
 *
 * @desc GPU Kernel to calculate the predicted position and velocity
 *       of all the particles.
 *
 * @note Working properly.
 */
__global__ void
k_predicted_pos_vel(double4 *d_r,   double4 *d_v,   double4 *d_a, double4 *d_a1,
                    double4 *d_p_r, double4 *d_p_v, double *d_t, double ITIME, int n)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    if (i < n)
    {
        double dt = ITIME - d_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        d_p_r[i].x = (dt3/6 * d_a1[i].x) + (dt2/2 * d_a[i].x) + (dt * d_v[i].x) + d_r[i].x;
        d_p_r[i].y = (dt3/6 * d_a1[i].y) + (dt2/2 * d_a[i].y) + (dt * d_v[i].y) + d_r[i].y;
        d_p_r[i].z = (dt3/6 * d_a1[i].z) + (dt2/2 * d_a[i].z) + (dt * d_v[i].z) + d_r[i].z;

        d_p_v[i].x = (dt2/2 * d_a1[i].x) + (dt  * d_a[i].x) + d_v[i].x;
        d_p_v[i].y = (dt2/2 * d_a1[i].y) + (dt  * d_a[i].y) + d_v[i].y;
        d_p_v[i].z = (dt2/2 * d_a1[i].z) + (dt  * d_a[i].z) + d_v[i].z;
    }
}
//__global__ void k_update_acc_jrk_single
//(double4 *new_a, double4 *new_j, double4 *r, double4 *v, float *m, int n, int current)
//{
//
//    double4 aa = {0.0, 0.0, 0.0, 0.0};
//    double4 jj = {0.0, 0.0, 0.0, 0.0};
//    double4 pos = r[current];
//    double4 vel = v[current];
//
//    int id = threadIdx.x + blockDim.x * blockIdx.x;
//
//    if (id < n)
//    {
//        if(id != current)
//        {
//            k_force_calculation(pos, vel, r[id], v[id], aa, jj, m[id]);
//        }
//
//        new_a[id] = aa;
//        new_j[id] = jj;
//    }
//}

//__global__ void k_update_2d(int *move,  double4 *new_acc, double4 *new_jrk,
//                            double4 *r, double4 *v,       float *m,
//                            int total,  int n)
//{
//
//    int xbid  = blockIdx.x;
//    int ybid  = blockIdx.y;
//    int tx  = threadIdx.x;
//    int gid = tx + blockDim.x * xbid;
//
//    int j_ini = (total * ybid)/NJBLOCK;
//    int j_end = (total * (ybid + 1))/NJBLOCK;
//
//    double4 pos = r[gid];
//    double4 vel = v[gid];
//    double4 j_acc  = {0.0,0.0,0.0,0.0};
//    double4 j_jrk = {0.0,0.0,0.0,0.0};
//
//    for (int j = j_ini; j < j_end; j+=BSIZE)
//    {
//        // Shared memory of BSIZE for j-particles
//        __shared__ double4 s_r[BSIZE];
//        __shared__ double4 s_v[BSIZE];
//        __syncthreads();
//
//        // Load of the r and v to shared memory of the j-particle
//        double4 *src_r = (double4 *)&r[j];
//        double4 *src_v = (double4 *)&v[j];
//        double4 *dst_r = (double4 *)s_r;
//        double4 *dst_v = (double4 *)s_v;
//
//        dst_r[tx]         = src_r[tx];
//        dst_r[BSIZE + tx] = src_r[BSIZE + tx];
//
//        dst_v[tx]         = src_v[tx];
//        dst_v[BSIZE + tx] = src_v[BSIZE + tx];
//
//        __syncthreads();
//
//        // If we need to work with an incomplete block
//        if (j_end - j < BSIZE)
//        {
//            #pragma unroll 4
//            for (int jj = 0; jj < j_end - j; jj++)
//            {
//                double4 pos_j = s_r[jj];
//                double4 vel_j = s_v[jj];
//                k_force_calculation(pos, vel, pos_j, vel_j, j_acc, j_jrk, m[jj]);
//            }
//        }
//        else
//        {
//            #pragma unroll 4
//            for (int jj = 0; jj < BSIZE; jj++)
//            {
//                double4 pos_j = s_r[jj];
//                double4 vel_j = s_v[jj];
//                k_force_calculation(pos, vel, pos_j, vel_j, j_acc, j_jrk, m[jj]);
//            }
//        }
//     }
//    new_acc[gid + n * ybid]  = j_acc;
//    new_jrk[gid + n * ybid] =  j_jrk;
//}
//
//__global__ void k_reduce(double4 *d_in, double4 *d_out, int n)
//{
//    extern __shared__ double4 sdata[];
//
//    unsigned int tid = threadIdx.x;
//    unsigned int i = blockIdx.x * (BSIZE * 2) + tid;
//    unsigned int gridSize = BSIZE * 2 * gridDim.x;
//    double4 sum = {0.0,0.0,0.0,0.0};
//
//    while (i < n)
//    {
//        sum = sum + d_in[i];
//        if (i + BSIZE < n)
//            sum = sum + d_in[i+BSIZE];
//        i += gridSize;
//    }
//
//    sdata[tid] = sum;
//    __syncthreads();
//
//    if(BSIZE >= 512){if(tid < 256){sdata[tid] += sdata[tid+256];} __syncthreads();}
//    if(BSIZE >= 256){if(tid < 128){sdata[tid] += sdata[tid+128];} __syncthreads();}
//    if(BSIZE >= 128){if(tid <  64){sdata[tid] += sdata[tid+ 64];} __syncthreads();}
//
//
//    if(tid < 32)
//    {
//        volatile double4* smem = sdata;
//        if (BSIZE >= 62) { smem[tid] += smem[tid + 32]; }
//        if (BSIZE >= 32) { smem[tid] += smem[tid + 16]; }
//        if (BSIZE >= 16) { smem[tid] += smem[tid +  8]; }
//        if (BSIZE >=  8) { smem[tid] += smem[tid +  4]; }
//        if (BSIZE >=  4) { smem[tid] += smem[tid +  2]; }
//        if (BSIZE >=  2) { smem[tid] += smem[tid +  1]; }
//    }
//
//    if(tid == 0) d_out[blockIdx.x] = sdata[0];
//}
//
