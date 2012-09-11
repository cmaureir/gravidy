#include "dynamics_gpu_kernels.cuh"

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
__global__
void
k_energy(double4 *r, double4 *v, float *ekin, float *epot, float *m, int n)
{
    float ekin_tmp = 0;
    float epot_tmp = 0;
    double r2, v2;
    double rx, ry, rz;
    double vx, vy, vz;

    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j;

    epot_tmp = 0;
    for (j = i+1; j < n; j++)
    {
        rx = r[j].x - r[i].x;
        ry = r[j].y - r[i].y;
        rz = r[j].z - r[i].z;

        r2 = rx*rx + ry*ry + rz*rz;// + E*E;

        epot_tmp -= (m[i] * m[j]) / sqrt(r2);
    }

    vx = v[i].x * v[i].x;
    vy = v[i].y * v[i].y;
    vz = v[i].z * v[i].z;

    v2 = vx + vy + vz;
    ekin_tmp = 0.5 * m[i] * v2;

    ekin[i] = ekin_tmp;
    epot[i] = epot_tmp;
}


/*
 * @fn k_init_acc_jek
 *
 * @param r    Position array, using only 3 dimensions, the last value of the double4 is empty.
 * @param v    Velocity array, using only 3 dimensions, the last value of the double4 is empty.
 * @param a    Acceleration array, using only 3 dimensions, the last value of the double4 is empty.
 * @param j    Jerk array, using only 3 dimensions, the last value of the double4 is empty.
 * @param m    Masses of the particles.
 * @param n    Number of particles.
 *
 * @brief
 *  CUDA Kernel to calculate the initial acceleration and jerk of the system.
 *
 */
__global__
void
k_init_acc_jerk(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int n)
{

    double3 tmp_r, tmp_v;
    double4 pos, vel;
    double4 tmp_a = {0.0, 0.0, 0.0, 0.0};
    double4 tmp_j = {0.0, 0.0, 0.0, 0.0};
    double f, f3, f5;
    float mj;

    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int k;

    if (i < n)
    {
        pos   = r[i];
        vel   = v[i];
        for (k = 0; k < n; k++)
        {
            mj = m[k];
            tmp_r.x = r[k].x - pos.x;
            tmp_r.y = r[k].y - pos.y;
            tmp_r.z = r[k].z - pos.z;

            tmp_v.x = v[k].x - vel.x;
            tmp_v.y = v[k].y - vel.y;
            tmp_v.z = v[k].z - vel.z;

            f = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + E*E;

            f3 = f * f * f;
            f5 = f3 * f * f;
            f3 = sqrt(f3);
            f5 = sqrt(f5);

            tmp_a.x += mj * tmp_r.x / f3;
            tmp_a.y += mj * tmp_r.y / f3;
            tmp_a.z += mj * tmp_r.z / f3;

            tmp_j.x += mj * (tmp_v.x/f3 + (3 * tmp_v.x * tmp_r.x * tmp_r.x)/f5);
            tmp_j.y += mj * (tmp_v.y/f3 + (3 * tmp_v.y * tmp_r.y * tmp_r.y)/f5);
            tmp_j.z += mj * (tmp_v.z/f3 + (3 * tmp_v.z * tmp_r.z * tmp_r.z)/f5);
        }
        a[i] = tmp_a;
        j[i] = tmp_j;
    }

}


__global__ void
k_init_acc_jerk_tile(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int n)
{

    extern __shared__ double4 sh[];
    double4 *s_r = (double4*)sh;
    double4 *s_v = (double4*)&s_r[blockDim.x];

    double3 tmp_r, tmp_v;
    double4 pos, vel;
    double4 tmp_a = {0.0, 0.0, 0.0, 0.0};
    double4 tmp_j = {0.0, 0.0, 0.0, 0.0};
    double f, f3, f5;
    float mj;

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int tx = threadIdx.x;
    int tile;

    if (id < n)
    {
        pos   = r[id];
        vel   = v[id];

        tile = 0;
        for (int i = 0; i < n; i += BSIZE)
        {
            int idx = tile * BSIZE + tx;
            
            s_r[tx]   = r[idx];
            s_v[tx]   = v[idx];
            mj = m[idx];
            __syncthreads();

            for (int k = 0; k < BSIZE; k++)
            {
                tmp_r.x = s_r[k].x - pos.x;
                tmp_r.y = s_r[k].y - pos.y;
                tmp_r.z = s_r[k].z - pos.z;

                tmp_v.x = s_v[k].x - vel.x;
                tmp_v.y = s_v[k].y - vel.y;
                tmp_v.z = s_v[k].z - vel.z;

                f = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + E*E;

                f3 = f * f * f;
                f5 = f3 * f * f;
                f3 = sqrt(f3);
                f5 = sqrt(f5);

                tmp_a.x += mj * tmp_r.x / f3;
                tmp_a.y += mj * tmp_r.y / f3;
                tmp_a.z += mj * tmp_r.z / f3;

                tmp_j.x += mj * (tmp_v.x/f3 + (3 * tmp_v.x * tmp_r.x * tmp_r.x)/f5);
                tmp_j.y += mj * (tmp_v.y/f3 + (3 * tmp_v.y * tmp_r.y * tmp_r.y)/f5);
                tmp_j.z += mj * (tmp_v.z/f3 + (3 * tmp_v.z * tmp_r.z * tmp_r.z)/f5);
            }
            __syncthreads();
            tile++;
        }

        a[id] = tmp_a;
        j[id] = tmp_j;
    }
}


/*
 * @fn k_update_acc_jerk()
 *
 * @brief
 *  GPU kernel that performs the force calculation for all the particles
 *
 * @todo
 *  Only perform the calculation of the particles with the t + dt = ITIME
 */
__global__ void
k_update_acc_jerk_simple(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int *move, int n, int total)
{

    double3 tmp_r;
    double3 tmp_v;
    double4 tmp_a = {0.0f,0.0f,0.0f,0.0f};
    double4 tmp_j = {0.0f,0.0f,0.0f,0.0f};
    double f, f3, f5;
    float mj;
    int  k;

    int id  = threadIdx.x + blockDim.x * blockIdx.x;
    int i = move[id];

    if(i != -1)
    { 
        for (k = 0; k < n; k++)
        {
            mj = m[k];
            if( i != k )
            {
                tmp_r.x = r[k].x - r[i].x;
                tmp_r.y = r[k].y - r[i].y;
                tmp_r.z = r[k].z - r[i].z;

                tmp_v.x = v[k].x - v[i].x;
                tmp_v.y = v[k].y - v[i].y;
                tmp_v.z = v[k].z - v[i].z;

                if(i == 0 || k == 0 || i == 1 || k == 1)
                {
                    #define NEW_E (1e-8)
                    f = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + NEW_E*NEW_E;
                }
                else
                {
                    f = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + E*E;
                }

                f3 = f * f * f;
                f5 = f3 * f * f;
                f3 = sqrt(f3);
                f5 = sqrt(f5);

                tmp_a.x += mj * tmp_r.x / f3;
                tmp_a.y += mj * tmp_r.y / f3;
                tmp_a.z += mj * tmp_r.z / f3;

                tmp_j.x += mj * (tmp_v.x/f3 + (3 * tmp_v.x * tmp_r.x * tmp_r.x)/f5);
                tmp_j.y += mj * (tmp_v.y/f3 + (3 * tmp_v.y * tmp_r.y * tmp_r.y)/f5);
                tmp_j.z += mj * (tmp_v.z/f3 + (3 * tmp_v.z * tmp_r.z * tmp_r.z)/f5);
            }
        }
        a[i] = tmp_a;
        j[i] = tmp_j;
    }

}

__global__ void
k_update_acc_jerk_tile(double4 *r, double4 *v, double4 *a, double4 *j, float *m, int *move, int n, int total)
{

    extern __shared__ double4 sh[];
    double4 *s_r = (double4*)sh;
    double4 *s_v = (double4*)&s_r[blockDim.x];

    double4 pos = {0.0f, 0.0f, 0.0f, 0.0f};
    double4 vel = {0.0f, 0.0f, 0.0f, 0.0f};
    double4 tmp_a = {0.0f, 0.0f, 0.0f, 0.0f};
    double4 tmp_j = {0.0f, 0.0f, 0.0f, 0.0f};

    double3 tmp_r,tmp_v;
    double f, f3, f5;
    float mj;

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

        s_r[tx] = r[idx];
        s_v[tx] = v[idx];
        mj      = m[idx];
        __syncthreads();

        if ( id != idx )
        {
            for (int k = 0; k < BSIZE; k++)
            {
                tmp_r.x = s_r[k].x - pos.x;
                tmp_r.y = s_r[k].y - pos.y;
                tmp_r.z = s_r[k].z - pos.z;

                tmp_v.x = s_v[k].x - vel.x;
                tmp_v.y = s_v[k].y - vel.y;
                tmp_v.z = s_v[k].z - vel.z;

                if(id == 0 || idx == 0 || id == 1 || idx == 1)
                {
                    #define NEW_E (1e-8)
                    f = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + NEW_E*NEW_E;
                }
                else
                {
                    f = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + E*E;
                }

                f3 = f * f * f;
                f5 = f3 * f * f;
                f3 = sqrt(f3);
                f5 = sqrt(f5);

                tmp_a.x += mj * tmp_r.x / f3;
                tmp_a.y += mj * tmp_r.y / f3;
                tmp_a.z += mj * tmp_r.z / f3;

                tmp_j.x += mj * (tmp_v.x/f3 + (3 * tmp_v.x * tmp_r.x * tmp_r.x)/f5);
                tmp_j.y += mj * (tmp_v.y/f3 + (3 * tmp_v.y * tmp_r.y * tmp_r.y)/f5);
                tmp_j.z += mj * (tmp_v.z/f3 + (3 * tmp_v.z * tmp_r.z * tmp_r.z)/f5);
            }
        }
        
        __syncthreads();
    }
    
    a[id] = tmp_a;
    j[id] = tmp_j;
}

__global__ void
k_correction_pos_vel(double4 *r,     double4 *v,     double4 *a, double4 *j,
                     double4 *old_a, double4 *old_j,
                     double4 *p_r,   double4 *p_v,   float  *t,  float  *dt, 
                     double ITIME,   int *move,      int total)

{

    // Acceleration derivates (2nd and 3rd)
    double ax0_2, ay0_2, az0_2;
    double ax0_3, ay0_3, az0_3;
    double ax1_2, ay1_2, az1_2;

    float dt1, dt2, dt3, dt4, dt5;
    double abs_a1, abs_j1, abs_j12, abs_a1_3, abs_a1_2, abs_a1_22;

    double tmp_dt;
    int i;

    int k = threadIdx.x + blockDim.x * blockIdx.x;
    
    if ( k < total)
    {
        i = move[k];

        t[i] = ITIME;
        dt1 = dt[i];
        dt2 = dt1 * dt1;
        dt3 = dt2 * dt1;
        dt4 = dt2 * dt2;
        dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        ax0_2 = (-6 * (old_a[i].x - a[i].x ) - dt1 * (4 * old_j[i].x + 2 * j[i].x) ) / dt2;
        ay0_2 = (-6 * (old_a[i].y - a[i].y ) - dt1 * (4 * old_j[i].y + 2 * j[i].y) ) / dt2;
        az0_2 = (-6 * (old_a[i].z - a[i].z ) - dt1 * (4 * old_j[i].z + 2 * j[i].z) ) / dt2;

        // Acceleration 3rd derivate
        ax0_3 = (12 * (old_a[i].x - a[i].x ) + 6 * dt1 * (old_j[i].x + j[i].x) ) / dt3;
        ay0_3 = (12 * (old_a[i].y - a[i].y ) + 6 * dt1 * (old_j[i].y + j[i].y) ) / dt3;
        az0_3 = (12 * (old_a[i].z - a[i].z ) + 6 * dt1 * (old_j[i].z + j[i].z) ) / dt3;

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
        ax1_2 = ax0_2 + dt1 * ax0_3;
        ay1_2 = ay0_2 + dt1 * ay0_3;
        az1_2 = az0_2 + dt1 * az0_3;

        // |a_{1,i}|
        abs_a1   = sqrt((a[i].x * a[i].x) + (a[i].y * a[i].y) + (a[i].z * a[i].z));
        // |j_{1,i}|
        abs_j1   = sqrt((j[i].x * j[i].x) + (j[i].y * j[i].y) + (j[i].z * j[i].z));
        // |j_{1,i}|^{2}
        abs_j12  = abs_j1 * abs_j1;
        // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
        abs_a1_3 = sqrt((ax0_3 * ax0_3) + (ay0_3 * ay0_3) + (az0_3 * az0_3));
        // |a_{1,i}^{(2)}|
        abs_a1_2 = sqrt((ax1_2 * ax1_2) + (ay1_2 * ay1_2) + (az1_2 * az1_2));
        // |a_{1,i}^{(2)}|^{2}
        abs_a1_22  = abs_a1_2 * abs_a1_2;

        tmp_dt = sqrt( ETA_N * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));

        /* Adjusting to block timesteps */
        if (tmp_dt < D_TIME_MIN)
        {
            tmp_dt = D_TIME_MIN;
        }
        else if (tmp_dt > D_TIME_MAX)
        {
            tmp_dt = D_TIME_MAX;
        }
        else
        {
            tmp_dt = powf(2,(int)((log(tmp_dt)/log(2.0))-1));
        }
        dt[i] = tmp_dt;
    }
}

__global__ void
k_update_acc_jerk_single(double4 current_pos, double4 current_vel, double4 *new_a, double4 *new_j, 
                         double4 *r,          double4 *v,          float *m,
                         int n,               int current)
{
    extern __shared__ double4 sh[];
    double4 *s_r = (double4*)sh;
    double4 *s_v = (double4*)&s_r[blockDim.x];

    double4 pos = current_pos;
    double4 vel = current_vel;

    double4 tmp_a = {0.0f, 0.0f, 0.0f, 0.0f};
    double4 tmp_j = {0.0f, 0.0f, 0.0f, 0.0f};

    double3 tmp_r,tmp_v;
    double f, f3, f5;
    float mj;

    int ii = threadIdx.x + blockDim.x * blockIdx.x;
    int tx = threadIdx.x;

    if (ii < n)
    {

        s_r[tx] = r[ii];
        s_v[tx] = v[ii];
        mj      = m[ii];
        __syncthreads();

        tmp_r.x = s_r[tx].x - pos.x;
        tmp_r.y = s_r[tx].y - pos.y;
        tmp_r.z = s_r[tx].z - pos.z;

        tmp_v.x = s_v[tx].x - vel.x;
        tmp_v.y = s_v[tx].y - vel.y;
        tmp_v.z = s_v[tx].z - vel.z;

        f = tmp_r.x*tmp_r.x + tmp_r.y*tmp_r.y + tmp_r.z*tmp_r.z + E*E;

        f3 = f * f * f;
        f5 = f3 * f * f;
        f3 = sqrt(f3);
        f5 = sqrt(f5);

        tmp_a.x = mj * tmp_r.x / f3;
        tmp_a.y = mj * tmp_r.y / f3;
        tmp_a.z = mj * tmp_r.z / f3;

        tmp_j.x = mj * (tmp_v.x/f3 + (3 * tmp_v.x * tmp_r.x * tmp_r.x)/f5);
        tmp_j.y = mj * (tmp_v.y/f3 + (3 * tmp_v.y * tmp_r.y * tmp_r.y)/f5);
        tmp_j.z = mj * (tmp_v.z/f3 + (3 * tmp_v.z * tmp_r.z * tmp_r.z)/f5);

        new_a[ii] = tmp_a;
        new_j[ii] = tmp_j;
    }
}
