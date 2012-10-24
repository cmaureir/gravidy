#include "dynamics_gpu_kernels.cuh"

inline __host__ __device__ double4 operator+(const double4 &a, const double4 &b)
{
    return make_double4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}

inline __host__ __device__ void operator+=(double4 &a, double4 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
}

inline __host__ __device__ void operator+=(volatile double4 &a, volatile double4 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;

}

__device__ void dev_gravity(double4 i_pos, double4 i_vel,
                            double4 j_pos, double4 j_vel,
                            double4 &acc,  double4 &jrk,
                            float j_mass)
{
    double dx =  j_pos.x - i_pos.x;
    double dy =  j_pos.y - i_pos.y;
    double dz =  j_pos.z - i_pos.z;
    double dvx = j_vel.x - i_vel.x;
    double dvy = j_vel.y - i_vel.y;
    double dvz = j_vel.z - i_vel.z;
    double r2 =  dx*dx + dy*dy + dz*dz;
    double rv =  dx*dvx + dy*dvy + dz*dvz;
    double rinv1 = rsqrtf(r2);
    double rinv2 = rinv1 * rinv1;
    double mrinv1 = j_mass * rinv1;
    double mrinv3 = mrinv1 * rinv2;
    rv *= -3.f * rinv2;
    acc.x += mrinv3 * dx;
    acc.y += mrinv3 * dy;
    acc.z += mrinv3 * dz;
    jrk.x += mrinv3 * (dvx + rv * dx);
    jrk.y += mrinv3 * (dvy + rv * dy);
    jrk.z += mrinv3 * (dvz + rv * dz);
}

__global__ void k_update_2d(int *move, double4 *new_acc, double4 *new_jrk,
                            double4 *r,       double4 *v,
                            float *m,         int total,  int n)
{

    // TO DO
    // Trabajar s√lo con el que recibimos
    int xbid  = blockIdx.x;
    int ybid  = blockIdx.y;
    int tx  = threadIdx.x;
    int gid = tx + blockDim.x * xbid;

    int j_ini = (total * ybid)/NJBLOCK;
    int j_end = (total * (ybid + 1))/NJBLOCK;

    double4 pos = r[gid];
    double4 vel = v[gid];
    double4 j_acc  = {0.0,0.0,0.0,0.0};
    double4 j_jrk = {0.0,0.0,0.0,0.0};

    for (int j = j_ini; j < j_end; j+=BSIZE)
    {
        // Shared memory of BSIZE for j-particles
        __shared__ double4 s_r[BSIZE];
        __shared__ double4 s_v[BSIZE];
        __syncthreads();

        // Load of the r and v to shared memory of the j-particle
        double4 *src_r = (double4 *)&r[j];
        double4 *src_v = (double4 *)&v[j];
        double4 *dst_r = (double4 *)s_r;
        double4 *dst_v = (double4 *)s_v;

        dst_r[tx]         = src_r[tx];
        dst_r[BSIZE + tx] = src_r[BSIZE + tx];

        dst_v[tx]         = src_v[tx];
        dst_v[BSIZE + tx] = src_v[BSIZE + tx];

        __syncthreads();

        // If we need to work with an incomplete block
        if (j_end - j < BSIZE)
        {
            #pragma unroll 4
            for (int jj = 0; jj < j_end - j; jj++)
            {
                double4 pos_j = s_r[jj];
                double4 vel_j = s_v[jj];
                dev_gravity(pos, vel, pos_j, vel_j, j_acc, j_jrk, m[jj]);
            }
        }
        else
        {
            #pragma unroll 4
            for (int jj = 0; jj < BSIZE; jj++)
            {
                double4 pos_j = s_r[jj];
                double4 vel_j = s_v[jj];
                dev_gravity(pos, vel, pos_j, vel_j, j_acc, j_jrk, m[jj]);
            }
        }
     }
    new_acc[gid + n * ybid]  = j_acc;
    new_jrk[gid + n * ybid] =  j_jrk;
}

__global__ void k_reduce_energy(float *d_in, float *d_out, int n)
{
    extern __shared__ float data[];

    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (BSIZE * 2) + tid;
    unsigned int gridSize = BSIZE * 2 * gridDim.x;
    float sum = 0.0f;

    while (i < n)
    {
        sum = sum + d_in[i];
        if (i + BSIZE < n)
            sum = sum + d_in[i+BSIZE];
        i += gridSize;
    }

    data[tid] = sum;
    __syncthreads();

    if(BSIZE >= 512){if(tid < 256){data[tid] += data[tid + 256];} __syncthreads();}
    if(BSIZE >= 256){if(tid < 128){data[tid] += data[tid + 128];} __syncthreads();}
    if(BSIZE >= 128){if(tid <  64){data[tid] += data[tid +  64];} __syncthreads();}

    if(tid < 32)
    {
        volatile float* smem = data;
        if (BSIZE >= 62) { smem[tid] += smem[tid + 32]; }
        if (BSIZE >= 32) { smem[tid] += smem[tid + 16]; }
        if (BSIZE >= 16) { smem[tid] += smem[tid +  8]; }
        if (BSIZE >=  8) { smem[tid] += smem[tid +  4]; }
        if (BSIZE >=  4) { smem[tid] += smem[tid +  2]; }
        if (BSIZE >=  2) { smem[tid] += smem[tid +  1]; }
    }

    if(tid == 0) d_out[blockIdx.x] = data[0];
}

__global__ void k_reduce(double4 *d_in, double4 *d_out, int n)
{
    extern __shared__ double4 sdata[];

    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (BSIZE * 2) + tid;
    unsigned int gridSize = BSIZE * 2 * gridDim.x;
    double4 sum = {0.0,0.0,0.0,0.0};

    while (i < n)
    {
        sum = sum + d_in[i];
        if (i + BSIZE < n)
            sum = sum + d_in[i+BSIZE];
        i += gridSize;
    }

    sdata[tid] = sum;
    __syncthreads();

    if(BSIZE >= 512){if(tid < 256){sdata[tid] += sdata[tid+256];} __syncthreads();}
    if(BSIZE >= 256){if(tid < 128){sdata[tid] += sdata[tid+128];} __syncthreads();}
    if(BSIZE >= 128){if(tid <  64){sdata[tid] += sdata[tid+ 64];} __syncthreads();}


    if(tid < 32)
    {
        volatile double4* smem = sdata;
        if (BSIZE >= 62) { smem[tid] += smem[tid + 32]; }
        if (BSIZE >= 32) { smem[tid] += smem[tid + 16]; }
        if (BSIZE >= 16) { smem[tid] += smem[tid +  8]; }
        if (BSIZE >=  8) { smem[tid] += smem[tid +  4]; }
        if (BSIZE >=  4) { smem[tid] += smem[tid +  2]; }
        if (BSIZE >=  2) { smem[tid] += smem[tid +  1]; }
    }

    if(tid == 0) d_out[blockIdx.x] = sdata[0];
}

/*
 * @fn k_energy
 *
 * @param r    Position array, using only 3 dimensions, the last value of the double4 is empty.
 * @param v    Velocity array, using only 3 dimensions, the last value of the double4 is empty.
 * @param ekin Kinetic energy array.
 * @param epot Potential energy array.
 * @param m    Masses of the particles.
 * @param n    Number of particles.
 *
 * @brief
 *  CUDA Kernel to calculate the potential and kinetic energy for each particle.
 *  The arrays contains the energy for each particle.
 *
 */
__global__ void
k_energy(double4 *r, double4 *v, float *ekin, float *epot, float *m, int n)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j;

    if (i < n)
    {
        float epot_tmp = 0.0f;
        for (j = i+1; j < n; j++)
        {
            double3 rr = {r[j].x - r[i].x, r[j].y - r[i].y, r[j].z - r[i].z};
            double r2 = rr.x*rr.x + rr.y*rr.y + rr.z*rr.z;

            epot_tmp -= (m[i] * m[j]) * rsqrt(r2);
        }

        double3 vv = {v[i].x * v[i].x, v[i].y * v[i].y, v[i].z * v[i].z};

        double v2 = vv.x + vv.y + vv.z;
        float ekin_tmp = 0.5 * m[i] * v2;

        ekin[i] = ekin_tmp;
        epot[i] = epot_tmp;
    }
}

__global__ void
k_init_acc_jrk(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int n)
{

    extern __shared__ double4 sh[];
    double4 *s_r = (double4*)sh;
    double4 *s_v = (double4*)&s_r[blockDim.x];

    double4 tmp_a = {0.0, 0.0, 0.0, 0.0};
    double4 tmp_j = {0.0, 0.0, 0.0, 0.0};
    float mj;

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int tx = threadIdx.x;

    if (id < n)
    {
        double4 pos   = r[id];
        double4 vel   = v[id];

        int tile = 0;
        for (int i = 0; i < n; i += BSIZE)
        {
            int idx = tile * BSIZE + tx;

            s_r[tx]   = r[idx];
            s_v[tx]   = v[idx];
            mj = m[idx];
            __syncthreads();

            for (int k = 0; k < BSIZE; k++)
            {
                double3 tmp_r = { s_r[k].x - pos.x, s_r[k].y - pos.y, s_r[k].z - pos.z};
                double3 tmp_v = { s_v[k].x - vel.x, s_v[k].y - vel.y, s_v[k].z - vel.z};

                double r2 = (double)(tmp_r.x*tmp_r.x + (tmp_r.y*tmp_r.y + (tmp_r.z*tmp_r.z + E2)));

                double rinv   = (double)rsqrt(r2);
                double r2inv  = (double)rinv  * rinv;
                double r3inv  = (double)r2inv * rinv;
                double r5inv  = (double)r2inv * r3inv;
                double mr3inv = (double)r3inv * mj;
                double mr5inv = (double)r5inv * mj;

                tmp_a.x += (double)(tmp_r.x * mr3inv);
                tmp_a.y += (double)(tmp_r.y * mr3inv);
                tmp_a.z += (double)(tmp_r.z * mr3inv);

                tmp_j.x += (double)(tmp_v.x * mr3inv + ((3 * (tmp_v.x * (tmp_r.x * tmp_r.x))) * mr5inv));
                tmp_j.y += (double)(tmp_v.y * mr3inv + ((3 * (tmp_v.y * (tmp_r.y * tmp_r.y))) * mr5inv));
                tmp_j.z += (double)(tmp_v.z * mr3inv + ((3 * (tmp_v.z * (tmp_r.z * tmp_r.z))) * mr5inv));
            }
            __syncthreads();
            tile++;
        }

        a[id] = tmp_a;
        j[id] = tmp_j;
    }
}


__global__ void k_update_acc_jrk(double4 *r, double4 *v, double4 *a, double4 *j,
                                 float *m,   int *move,  int n,      int total)
{

    extern __shared__ double4 sh[];
    double4 *s_r = (double4*)sh;
    double4 *s_v = (double4*)&s_r[blockDim.x];

    double4 pos = {0.0, 0.0, 0.0, 0.0};
    double4 vel = {0.0, 0.0, 0.0, 0.0};
    double4 tmp_a = {0.0, 0.0, 0.0, 0.0};
    double4 tmp_j = {0.0, 0.0, 0.0, 0.0};

    int ii = threadIdx.x + blockDim.x * blockIdx.x;
    int tx = threadIdx.x;
    int id = move[ii];

    if (id != -1)
    {
        pos   = r[id];
        vel   = v[id];
    }

    for (int i = 0; i < n; i += BSIZE)
    {
        int idx = i + tx;

        s_r[tx]  = r[idx];
        s_v[tx]  = v[idx];
        float mj = m[idx];
        __syncthreads();

        for (int k = 0; k < BSIZE; k++)
        {
            if(id != idx)
            {
                double3 tmp_r = { s_r[k].x - pos.x, s_r[k].y - pos.y, s_r[k].z - pos.z};
                double3 tmp_v = { s_v[k].x - vel.x, s_v[k].y - vel.y, s_v[k].z - vel.z};

                double r2 = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + E2;

                double rinv = rsqrt(r2);
                double r2inv = rinv  * rinv;
                double r3inv = r2inv * rinv;
                double r5inv = r2inv * r3inv;
                double mr3inv = r3inv * mj;
                double mr5inv = r5inv * mj;

                tmp_a.x += tmp_r.x * mr3inv;
                tmp_a.y += tmp_r.y * mr3inv;
                tmp_a.z += tmp_r.z * mr3inv;

                tmp_j.x += tmp_v.x * mr3inv + (3 * tmp_v.x * tmp_r.x * tmp_r.x) * mr5inv;
                tmp_j.y += tmp_v.y * mr3inv + (3 * tmp_v.y * tmp_r.y * tmp_r.y) * mr5inv;
                tmp_j.z += tmp_v.z * mr3inv + (3 * tmp_v.z * tmp_r.z * tmp_r.z) * mr5inv;
            }
        }
        __syncthreads();
    }

    if (id != -1)
    {
        a[id] = tmp_a;
        j[id] = tmp_j;
    }
}

__global__ void
k_correction_pos_vel(double4 *r,     double4 *v,     double4 *a,   double4 *j,
                     double4 *old_a, double4 *old_j, double4 *p_r, double4 *p_v,
                     float  *t,      float  *dt,     float ITIME, int *move,
                     int total)
{

    int k = threadIdx.x + blockDim.x * blockIdx.x;

    if ( k < total)
    {
        int i = move[k];

        float dt1 = ITIME - t[i];
        t[i] = ITIME;
        float dt2 = dt1 * dt1;
        float dt3 = dt2 * dt1;
        float dt4 = dt2 * dt2;
        float dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        double ax0_2 = (-6 * (old_a[i].x - a[i].x ) - dt1 * (4 * old_j[i].x + 2 * j[i].x) ) / dt2;
        double ay0_2 = (-6 * (old_a[i].y - a[i].y ) - dt1 * (4 * old_j[i].y + 2 * j[i].y) ) / dt2;
        double az0_2 = (-6 * (old_a[i].z - a[i].z ) - dt1 * (4 * old_j[i].z + 2 * j[i].z) ) / dt2;

        // Acceleration 3rd derivate
        double ax0_3 = (12 * (old_a[i].x - a[i].x ) + 6 * dt1 * (old_j[i].x + j[i].x) ) / dt3;
        double ay0_3 = (12 * (old_a[i].y - a[i].y ) + 6 * dt1 * (old_j[i].y + j[i].y) ) / dt3;
        double az0_3 = (12 * (old_a[i].z - a[i].z ) + 6 * dt1 * (old_j[i].z + j[i].z) ) / dt3;

        // Correcting position
        r[i].x = p_r[i].x + (dt4/24)*ax0_2 + (dt5/120)*ax0_3;
        r[i].y = p_r[i].y + (dt4/24)*ay0_2 + (dt5/120)*ay0_3;
        r[i].z = p_r[i].z + (dt4/24)*az0_2 + (dt5/120)*az0_3;

        // Correcting velocity
        v[i].x = p_v[i].x + (dt3/6)*ax0_2 + (dt4/24)*ax0_3;
        v[i].y = p_v[i].y + (dt3/6)*ay0_2 + (dt4/24)*ay0_3;
        v[i].z = p_v[i].z + (dt3/6)*az0_2 + (dt4/24)*az0_3;

        // Timestep update

        // Calculating a_{1,i}^{(2)} = a_{0,i}^{(2)} + dt * a_{0,i}^{(3)}
        double ax1_2 = ax0_2 + dt1 * ax0_3;
        double ay1_2 = ay0_2 + dt1 * ay0_3;
        double az1_2 = az0_2 + dt1 * az0_3;

        // |a_{1,i}|
        double abs_a1   = sqrt((a[i].x * a[i].x) + (a[i].y * a[i].y) + (a[i].z * a[i].z));
        // |j_{1,i}|
        double abs_j1   = sqrt((j[i].x * j[i].x) + (j[i].y * j[i].y) + (j[i].z * j[i].z));
        // |j_{1,i}|^{2}
        double abs_j12  = abs_j1 * abs_j1;
        // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
        double abs_a1_3 = sqrt((ax0_3 * ax0_3) + (ay0_3 * ay0_3) + (az0_3 * az0_3));
        // |a_{1,i}^{(2)}|
        double abs_a1_2 = sqrt((ax1_2 * ax1_2) + (ay1_2 * ay1_2) + (az1_2 * az1_2));
        // |a_{1,i}^{(2)}|^{2}
        double abs_a1_22  = abs_a1_2 * abs_a1_2;

        double tmp_dt = sqrt(ETA_N * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));

        /* Adjusting to block timesteps */
        if (tmp_dt < dt[i])
            tmp_dt = dt[i]/2;
        else if (tmp_dt > 2 * dt[i])
            tmp_dt = 2 * dt[i];
        else
            tmp_dt = dt[i];

        if (tmp_dt < D_TIME_MIN)
            tmp_dt = D_TIME_MIN;
        else if (tmp_dt > D_TIME_MAX)
            tmp_dt = D_TIME_MAX;

        dt[i] = tmp_dt;
    }
}

__global__ void
k_update_acc_jrk_single(double4 c_pos, double4 c_vel, double4 *new_a, double4 *new_j,
                         double4 *r,          double4 *v,          float *m,
                         int n,               int current)
{

    double4 tmp_a = {0.0, 0.0, 0.0, 0.0};
    double4 tmp_j = {0.0, 0.0, 0.0, 0.0};

    int ii = threadIdx.x + blockDim.x * blockIdx.x;

    if (ii < n)
    {
        new_a[ii] = tmp_a;
        new_j[ii] = tmp_j;
        if(ii != current)
        {
            float mj = m[ii];
            double3 tmp_r = { r[ii].x - c_pos.x, r[ii].y - c_pos.y, r[ii].z - c_pos.z};
            double3 tmp_v = { v[ii].x - c_vel.x, v[ii].y - c_vel.y, v[ii].z - c_vel.z};

            double r2 = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + E2;

            double rinv = rsqrt(r2);
            double r2inv = rinv  * rinv;
            double r3inv = r2inv * rinv;
            double r5inv = r2inv * r3inv;
            double mr3inv = r3inv * mj;
            double mr5inv = r5inv * mj;

            tmp_a.x = tmp_r.x * mr3inv;
            tmp_a.y = tmp_r.y * mr3inv;
            tmp_a.z = tmp_r.z * mr3inv;

            tmp_j.x = tmp_v.x * mr3inv + (3 * tmp_v.x * tmp_r.x * tmp_r.x) * mr5inv;
            tmp_j.y = tmp_v.y * mr3inv + (3 * tmp_v.y * tmp_r.y * tmp_r.y) * mr5inv;
            tmp_j.z = tmp_v.z * mr3inv + (3 * tmp_v.z * tmp_r.z * tmp_r.z) * mr5inv;

            new_a[ii] = tmp_a;
            new_j[ii] = tmp_j;
        }
    }
}

__global__ void
k_predicted_pos_vel(double4 *d_r,   double4 *d_v,   double4 *d_a, double4 *d_j,
                    double4 *d_p_r, double4 *d_p_v, float *d_t, float ITIME, int n)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    if (i < n)
    {
        float dt = ITIME - d_t[i];
        float dt2 = (dt  * dt)/2;
        float dt3 = (dt2 * dt)/6;

        d_p_r[i].x = dt3 * d_j[i].x + dt2 * d_a[i].x + dt * d_v[i].x + d_r[i].x;
        d_p_r[i].y = dt3 * d_j[i].y + dt2 * d_a[i].y + dt * d_v[i].y + d_r[i].y;
        d_p_r[i].z = dt3 * d_j[i].z + dt2 * d_a[i].z + dt * d_v[i].z + d_r[i].z;

        d_p_v[i].x = dt2 * d_j[i].x + dt * d_a[i].x + d_v[i].x;
        d_p_v[i].y = dt2 * d_j[i].y + dt * d_a[i].y + d_v[i].y;
        d_p_v[i].z = dt2 * d_j[i].z + dt * d_a[i].z + d_v[i].z;
    }
}

__global__ void
k_update_acc_jrk_simple(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int *move, int n, int total)
{

    double4 tmp_a = {0.0, 0.0, 0.0, 0.0};
    double4 tmp_j = {0.0, 0.0, 0.0, 0.0};

    int id = threadIdx.x + blockDim.x * blockIdx.x;

    if (id < n)
    {
        double4 pos   = r[id];
        double4 vel   = v[id];
        float   mj    = m[id];

        for (int i = 0; i < n; i++)
        {
            if(i == id) continue;

            double3 tmp_r = { r[i].x - pos.x, r[i].y - pos.y, r[i].z - pos.z};
            double3 tmp_v = { v[i].x - vel.x, v[i].y - vel.y, v[i].z - vel.z};

            double r2     = (double)(tmp_r.x*tmp_r.x + (tmp_r.y*tmp_r.y + (tmp_r.z*tmp_r.z + E2)));

            double rinv   = (double)rsqrt(r2);
            double r2inv  = (double)rinv  * rinv;
            double r3inv  = (double)r2inv * rinv;
            double r5inv  = (double)r2inv * r3inv;
            double mr3inv = (double)r3inv * mj;
            double mr5inv = (double)r5inv * mj;

            tmp_a.x += (double)(tmp_r.x * mr3inv);
            tmp_a.y += (double)(tmp_r.y * mr3inv);
            tmp_a.z += (double)(tmp_r.z * mr3inv);

            tmp_j.x += (double)(tmp_v.x * mr3inv + ((3 * (tmp_v.x * (tmp_r.x * tmp_r.x))) * mr5inv));
            tmp_j.y += (double)(tmp_v.y * mr3inv + ((3 * (tmp_v.y * (tmp_r.y * tmp_r.y))) * mr5inv));
            tmp_j.z += (double)(tmp_v.z * mr3inv + ((3 * (tmp_v.z * (tmp_r.z * tmp_r.z))) * mr5inv));
        }

        a[id] = tmp_a;
        j[id] = tmp_j;
    }
}
