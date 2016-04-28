#include "Hermite4GPU.cuh"

Hermite4GPU::~Hermite4GPU()
{
    free_arrays_device();
}

void Hermite4GPU::alloc_arrays_device()
{

    size_t i1_size = ns->n * sizeof(int);
    size_t d1_size = ns->n * sizeof(double);
    size_t d4_size = ns->n * sizeof(double4);
    size_t ff_size = ns->n * sizeof(Forces);
    size_t pp_size = ns->n * sizeof(Predictor);

    for(int g = 0; g < gpus; g++)
    {
        // Setting GPU
        CSC(cudaSetDevice(g));

        CSC(cudaMalloc((void**)&ns->d_r[g], d4_size));
        CSC(cudaMalloc((void**)&ns->d_v[g], d4_size));
        CSC(cudaMalloc((void**)&ns->d_f[g], ff_size));
        CSC(cudaMalloc((void**)&ns->d_p[g], pp_size));
        CSC(cudaMalloc((void**)&ns->d_t[g], d1_size));
        CSC(cudaMalloc((void**)&ns->d_i[g], pp_size));
        CSC(cudaMalloc((void**)&ns->d_dt[g], d1_size));
        CSC(cudaMalloc((void**)&ns->d_ekin[g], d1_size));
        CSC(cudaMalloc((void**)&ns->d_epot[g], d1_size));
        CSC(cudaMalloc((void**)&ns->d_move[g], i1_size));
        CSC(cudaMalloc((void**)&ns->d_fout[g], ff_size * NJBLOCK));
        CSC(cudaMalloc((void**)&ns->d_fout_tmp[g], ff_size * NJBLOCK));

        CSC(cudaMemset(ns->d_r[g], 0, d4_size));
        CSC(cudaMemset(ns->d_v[g], 0, d4_size));
        CSC(cudaMemset(ns->d_f[g], 0, ff_size));
        CSC(cudaMemset(ns->d_p[g], 0, pp_size));
        CSC(cudaMemset(ns->d_t[g], 0, d1_size));
        CSC(cudaMemset(ns->d_i[g], 0, pp_size));
        CSC(cudaMemset(ns->d_dt[g], 0, d1_size));
        CSC(cudaMemset(ns->d_ekin[g], 0, d1_size));
        CSC(cudaMemset(ns->d_epot[g], 0, d1_size));
        CSC(cudaMemset(ns->d_move[g], 0, i1_size));
        CSC(cudaMemset(ns->d_fout[g], 0, ff_size * NJBLOCK));
        CSC(cudaMemset(ns->d_fout_tmp[g], 0, ff_size * NJBLOCK));

        ns->h_fout_gpu[g] = new Forces[ns->n*NJBLOCK];
    }

    // Extra CPU array
    ns->h_fout_tmp = new Forces[ns->n*NJBLOCK];

}

void Hermite4GPU::free_arrays_device()
{

    for(int g = 0; g < gpus; g++)
    {
        // Setting GPU
        CSC(cudaSetDevice(g));

        CSC(cudaFree(ns->d_r[g]));
        CSC(cudaFree(ns->d_v[g]));
        CSC(cudaFree(ns->d_f[g]));
        CSC(cudaFree(ns->d_p[g]));
        CSC(cudaFree(ns->d_t[g]));
        CSC(cudaFree(ns->d_i[g]));
        CSC(cudaFree(ns->d_dt[g]));
        CSC(cudaFree(ns->d_ekin[g]));
        CSC(cudaFree(ns->d_epot[g]));
        CSC(cudaFree(ns->d_move[g]));
        CSC(cudaFree(ns->d_fout[g]));
        CSC(cudaFree(ns->d_fout_tmp[g]));
        delete ns->h_fout_gpu[g];
    }

    delete ns->h_fout_tmp;
    //delete ns->h_fout_gpu;
}

/** Not implemented using GPU */
void Hermite4GPU::predicted_pos_vel(double ITIME)
{

    ns->gtime.prediction_ini = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < ns->n; i++)
    {
        double dt  = ITIME - ns->h_t[i];
        double dt2 = 0.5*(dt  * dt);
        double dt3 = 0.166666666666666*(dt * dt * dt);

        Forces ff = ns->h_f[i];
        double4 rr = ns->h_r[i];
        double4 vv = ns->h_v[i];

        ns->h_p[i].r[0] = (dt3 * ff.a1[0]) + (dt2 * ff.a[0]) + (dt * vv.x) + rr.x;
        ns->h_p[i].r[1] = (dt3 * ff.a1[1]) + (dt2 * ff.a[1]) + (dt * vv.y) + rr.y;
        ns->h_p[i].r[2] = (dt3 * ff.a1[2]) + (dt2 * ff.a[2]) + (dt * vv.z) + rr.z;

        ns->h_p[i].v[0] = (dt2 * ff.a1[0]) + (dt * ff.a[0]) + vv.x;
        ns->h_p[i].v[1] = (dt2 * ff.a1[1]) + (dt * ff.a[1]) + vv.y;
        ns->h_p[i].v[2] = (dt2 * ff.a1[2]) + (dt * ff.a[2]) + vv.z;

        ns->h_p[i].m = rr.w;

    }
    ns->gtime.prediction_end += omp_get_wtime() - ns->gtime.prediction_ini;
}

/** Not implemented using GPU */
void Hermite4GPU::correction_pos_vel(double ITIME, int nact)
{
    ns->gtime.correction_ini = omp_get_wtime();
    int i;
    #pragma omp parallel for private(i)
    for (int k = 0; k < nact; k++)
    {
        i = ns->h_move[k];

        Forces ff = ns->h_f[i];
        Forces oo = ns->h_old[i];
        Predictor pp = ns->h_p[i];

        double dt1 = ns->h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        double dt2inv = 1.0/dt2;
        double dt3inv = 1.0/dt3;

        double dt3_6 = 0.166666666666666*dt3;
        double dt4_24 = 0.041666666666666*dt4;
        double dt5_120 = 0.008333333333333*dt5;

        // Acceleration 2nd derivate
        ns->h_a2[i].x = (-6 * (oo.a[0] - ff.a[0] ) - dt1 * (4 * oo.a1[0] + 2 * ff.a1[0]) ) * dt2inv;
        ns->h_a2[i].y = (-6 * (oo.a[1] - ff.a[1] ) - dt1 * (4 * oo.a1[1] + 2 * ff.a1[1]) ) * dt2inv;
        ns->h_a2[i].z = (-6 * (oo.a[2] - ff.a[2] ) - dt1 * (4 * oo.a1[2] + 2 * ff.a1[2]) ) * dt2inv;

        // Acceleration 3rd derivate
        ns->h_a3[i].x = (12 * (oo.a[0] - ff.a[0] ) + 6 * dt1 * (oo.a1[0] + ff.a1[0]) ) * dt3inv;
        ns->h_a3[i].y = (12 * (oo.a[1] - ff.a[1] ) + 6 * dt1 * (oo.a1[1] + ff.a1[1]) ) * dt3inv;
        ns->h_a3[i].z = (12 * (oo.a[2] - ff.a[2] ) + 6 * dt1 * (oo.a1[2] + ff.a1[2]) ) * dt3inv;


        // Correcting position
        ns->h_r[i].x = pp.r[0] + (dt4_24)*ns->h_a2[i].x + (dt5_120)*ns->h_a3[i].x;
        ns->h_r[i].y = pp.r[1] + (dt4_24)*ns->h_a2[i].y + (dt5_120)*ns->h_a3[i].y;
        ns->h_r[i].z = pp.r[2] + (dt4_24)*ns->h_a2[i].z + (dt5_120)*ns->h_a3[i].z;

        // Correcting velocity
        ns->h_v[i].x = pp.v[0] + (dt3_6)*ns->h_a2[i].x + (dt4_24)*ns->h_a3[i].x;
        ns->h_v[i].y = pp.v[1] + (dt3_6)*ns->h_a2[i].y + (dt4_24)*ns->h_a3[i].y;
        ns->h_v[i].z = pp.v[2] + (dt3_6)*ns->h_a2[i].z + (dt4_24)*ns->h_a3[i].z;

        ns->h_t[i] = ITIME;

        double normal_dt  = nu->get_timestep_normal(i, ns->eta);
        ns->h_dt[i] = nu->normalize_dt(normal_dt, ns->h_dt[i], ns->h_t[i], i);

    }
    ns->gtime.correction_end += omp_get_wtime() - ns->gtime.correction_ini;
}

void Hermite4GPU::init_acc_jrk()
{

    size_t pp_size = ns->n * sizeof(Predictor);

    // Copying arrays to device
    for(int g = 0; g < gpus; g++)
    {
        CSC(cudaSetDevice(g));

        // All this information from the predictors is needed by each device
        CSC(cudaMemcpy(ns->d_p[g], ns->h_p, pp_size, cudaMemcpyHostToDevice));
    }

    // Executing kernels
    for(int g = 0; g < gpus; g++)
    {
        CSC(cudaSetDevice(g));

        nthreads = BSIZE;
        nblocks = std::ceil(n_part[g]/(float)nthreads);

        k_init_acc_jrk <<< nblocks, nthreads, smem >>> (ns->d_p[g],
                                                        ns->d_f[g],
                                                        ns->n,
                                                        ns->e2,
                                                        g,
                                                        n_part[g]);
        get_kernel_error();
    }

    for(int g = 0; g < gpus; g++)
    {
        CSC(cudaSetDevice(g));

        size_t chunk = n_part[g]*sizeof(Forces);
        size_t slice = g*n_part[g-1];

        CSC(cudaMemcpy(&ns->h_f[slice], ns->d_f[g], chunk, cudaMemcpyDeviceToHost));
    }
}

void Hermite4GPU::update_acc_jrk(int nact)
{
    // Timer begin
    ns->gtime.update_ini = omp_get_wtime();

    for(int g = 0; g < gpus; g++)
    {
        if (n_part[g] > 0)
        {
            size_t pp_size = n_part[g] * sizeof(Predictor);
            int shift = g*n_part[g-1];

            CSC(cudaSetDevice(g));
            // Copying to the device the predicted r and v
            CSC(cudaMemcpy(ns->d_p[g], ns->h_p + shift, pp_size, cudaMemcpyHostToDevice));
        }
    }

    // Fill the h_i Predictor array with the particles that we need to move
    #pragma omp parallel for
    for (int i = 0; i < nact; i++)
    {
        ns->h_i[i] = ns->h_p[ns->h_move[i]];
    }

    for(int g = 0; g < gpus; g++)
    {
        if (n_part[g] > 0)
        {
            CSC(cudaSetDevice(g));
            // Copy to the GPU (d_i) the preddictor host array (h_i)
            size_t chunk = nact * sizeof(Predictor);
            CSC(cudaMemcpy(ns->d_i[g], ns->h_i, chunk, cudaMemcpyHostToDevice));
        }
    }

    ns->gtime.grav_ini = omp_get_wtime();
    for(int g = 0; g < gpus; g++)
    {
        if (n_part[g] > 0)
        {
            CSC(cudaSetDevice(g));
            // Blocks, threads and shared memory configuration
            int  nact_blocks = 1 + (nact-1)/BSIZE;
            dim3 nblocks(nact_blocks, NJBLOCK, 1);
            dim3 nthreads(BSIZE, 1, 1);

            // Kernel to update the forces for the particles in d_i
            k_update <<< nblocks, nthreads, smem >>> (ns->d_i[g],
                                                      ns->d_p[g], // partial
                                                      ns->d_fout[g],
                                                      n_part[g], // former N
                                                      nact,
                                                      ns->e2);
        }
    }

    ns->gtime.grav_end += omp_get_wtime() - ns->gtime.grav_ini;
    get_kernel_error();

    ns->gtime.reduce_ini = omp_get_wtime();
    for(int g = 0; g < gpus; g++)
    {
        if (n_part[g] > 0)
        {
            CSC(cudaSetDevice(g));
            // Blocks, threads and shared memory configuration for the reduction.
            dim3 rgrid   (nact,   1, 1);
            dim3 rthreads(NJBLOCK, 1, 1);

            // Kernel to reduce que temp array with the forces
            reduce <<< rgrid, rthreads, smem_reduce >>>(ns->d_fout[g],
                                                        ns->d_fout_tmp[g]);
        }
    }

    ns->gtime.reduce_end += omp_get_wtime() - ns->gtime.reduce_ini;
    get_kernel_error();

    for(int g = 0; g < gpus; g++)
    {
        if (n_part[g] > 0)
        {
            CSC(cudaSetDevice(g));
            size_t chunk = nact*sizeof(Forces);

            // Copy from the GPU the new forces for the d_i particles.
            CSC(cudaMemcpy(ns->h_fout_gpu[g], ns->d_fout_tmp[g], chunk,
                              cudaMemcpyDeviceToHost));
        }
    }



    // Update forces in the host
    ns->gtime.reduce_forces_ini = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < nact; i++)
    {
        int id = ns->h_move[i];
        ns->h_f[id].a[0] = 0.0;
        ns->h_f[id].a[1] = 0.0;
        ns->h_f[id].a[2] = 0.0;
        ns->h_f[id].a1[0] = 0.0;
        ns->h_f[id].a1[1] = 0.0;
        ns->h_f[id].a1[2] = 0.0;

        for(int g = 0; g < gpus; g++)
        {
            if (n_part[g] > 0)
            {
                ns->h_f[id] += ns->h_fout_gpu[g][i];
            }
        }
    }
    ns->gtime.reduce_forces_end += omp_get_wtime() - ns->gtime.reduce_forces_ini;

    // Timer end
    ns->gtime.update_end += (omp_get_wtime() - ns->gtime.update_ini);
}

double Hermite4GPU::get_energy_gpu()
{

    double time_energy_ini = omp_get_wtime();

    for(int g = 0; g < gpus; g++)
    {
        CSC(cudaSetDevice(g));

        size_t d4_size = ns->n * sizeof(double4);

        CSC(cudaMemcpy(ns->d_r[g], ns->h_r, d4_size, cudaMemcpyHostToDevice));
        CSC(cudaMemcpy(ns->d_v[g], ns->h_v, d4_size, cudaMemcpyHostToDevice));
    }

    for(int g = 0; g < gpus; g++)
    {
        CSC(cudaSetDevice(g));

        int nthreads = BSIZE;
        int nblocks = std::ceil(n_part[g]/(float)nthreads);
        k_energy <<< nblocks, nthreads >>> (ns->d_r[g],
                                            ns->d_v[g],
                                            ns->d_ekin[g],
                                            ns->d_epot[g],
                                            ns->n,
                                            n_part[g],
                                            g);
    }

    for(int g = 0; g < gpus; g++)
    {
        CSC(cudaSetDevice(g));

        size_t chunk = n_part[g]*sizeof(double);
        size_t slice = g*n_part[g-1];

        CSC(cudaMemcpy(&ns->h_ekin[slice], ns->d_ekin[g], chunk, cudaMemcpyDeviceToHost));
        CSC(cudaMemcpy(&ns->h_epot[slice], ns->d_epot[g], chunk, cudaMemcpyDeviceToHost));
    }

    // Reduction on CPU
    ns->en.kinetic = 0.0;
    ns->en.potential = 0.0;

    for (int i = 0; i < ns->n; i++)
    {
        ns->en.kinetic   += ns->h_ekin[i];
        ns->en.potential += ns->h_epot[i];
    }
    double time_energy_end = omp_get_wtime() - time_energy_ini;
    return ns->en.kinetic + ns->en.potential;
}

/*
 * @fn k_init_acc_jr
 *
 *
 * @desc GPU Kernel which calculates the initial acceleration and jerk
 * of all the particles of the system.
 *
 */
__global__ void k_init_acc_jrk (Predictor *p,
                                Forces *f,
                                int n,
                                double e2,
                                int dev,
                                int dev_size)
{

    extern __shared__ Predictor sh[];

    Forces ff;
    ff.a[0]  = 0.0;
    ff.a[1]  = 0.0;
    ff.a[2]  = 0.0;
    ff.a1[0] = 0.0;
    ff.a1[1] = 0.0;
    ff.a1[2] = 0.0;

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int tx = threadIdx.x;

    if (id < dev_size)
    {
      Predictor pred = p[id+(dev*dev_size)];
      //Predictor pred = p[id];
      int tile = 0;
      for (int i = 0; i < n; i += BSIZE)
      {
          int idx = tile * BSIZE + tx;
          sh[tx]   = p[idx];
          __syncthreads();
          for (int k = 0; k < BSIZE; k++)
          {
              k_force_calculation(pred, sh[k], ff, e2);
          }
          __syncthreads();
          tile++;
      }
      f[id] = ff;
    }
}

__device__ void k_force_calculation(Predictor i_p,
                                    Predictor j_p,
                                    Forces &f,
                                    double e2)
{
    double rx = j_p.r[0] - i_p.r[0];
    double ry = j_p.r[1] - i_p.r[1];
    double rz = j_p.r[2] - i_p.r[2];

    double vx = j_p.v[0] - i_p.v[0];
    double vy = j_p.v[1] - i_p.v[1];
    double vz = j_p.v[2] - i_p.v[2];

    double r2     = rx*rx + ry*ry + rz*rz + e2;
    double rinv   = rsqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * j_p.m;
    double mr5inv = r5inv * j_p.m;

    double rv = rx*vx + ry*vy + rz*vz;

    f.a[0] += (rx * mr3inv);
    f.a[1] += (ry * mr3inv);
    f.a[2] += (rz * mr3inv);

    f.a1[0] += (vx * mr3inv - (3 * rv) * rx * mr5inv);
    f.a1[1] += (vy * mr3inv - (3 * rv) * ry * mr5inv);
    f.a1[2] += (vz * mr3inv - (3 * rv) * rz * mr5inv);
}

/*
 * @fn k_update()
 *
 * @brief Gravitational interaction kernel.
 */
__global__ void k_update(Predictor *i_p,
                         Predictor *j_p,
                         Forces *fout,
                         int n,
                         int total,
                         double e2)
{
    int ibid = blockIdx.x;
    int jbid = blockIdx.y;
    int tid  = threadIdx.x;
    int iaddr  = tid + blockDim.x * ibid;
    int jstart = (n * (jbid  )) / NJBLOCK;
    int jend   = (n * (jbid+1)) / NJBLOCK;

    Predictor ip = i_p[iaddr];
    Forces fo;
    fo.a[0] = 0.0;
    fo.a[1] = 0.0;
    fo.a[2] = 0.0;
    fo.a1[0] = 0.0;
    fo.a1[1] = 0.0;
    fo.a1[2] = 0.0;

        for(int j=jstart; j<jend; j+=BSIZE)
        {
            __shared__ Predictor jpshare[BSIZE];
            __syncthreads();
            Predictor *src = (Predictor *)&j_p[j];
            Predictor *dst = (Predictor *)jpshare;
            dst[      tid] = src[      tid];
            dst[BSIZE+tid] = src[BSIZE+tid];
            __syncthreads();

            // If the total amount of particles is not a multiple of BSIZE
            if(jend-j < BSIZE)
            {
                #pragma unroll 4
                for(int jj=0; jj<jend-j; jj++)
                {
                    Predictor jp = jpshare[jj];
                    k_force_calculation(ip, jp, fo, e2);
                }
            }
            else
            {
                #pragma unroll 4
                for(int jj=0; jj<BSIZE; jj++)
                {
                    Predictor jp = jpshare[jj];
                    k_force_calculation(ip, jp, fo, e2);
                }
            }
        }
        fout[iaddr*NJBLOCK + jbid] = fo;
}

/*
 * @fn reduce()
 *
 * @brief Forces reduction kernel
 */
__global__ void reduce(Forces *in,
                       Forces *out)
{
    extern __shared__ Forces sdata[];

    const int xid   = threadIdx.x;
    const int bid   = blockIdx.x;
    const int iaddr = xid + blockDim.x * bid;

    sdata[xid] = in[iaddr];
    __syncthreads();

    if(xid < 8) sdata[xid] += sdata[xid + 8];
    if(xid < 4) sdata[xid] += sdata[xid + 4];
    if(xid < 2) sdata[xid] += sdata[xid + 2];
    if(xid < 1) sdata[xid] += sdata[xid + 1];

    if(xid == 0){
        out[bid] = sdata[0];
    }
}

__global__ void k_energy(double4 *r,
                         double4 *v,
                         double *ekin,
                         double *epot,
                         int n,
                         int dev_size,
                         int dev)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j;
    double ekin_tmp = 0.0;
    int id = i+dev*dev_size;

    if (i < dev_size)
    {
        double epot_tmp = 0.0;
        double4 ri = r[id];
        double4 vi = v[id];
        for (j = id+1; j < n; j++)
        {
            double rx = r[j].x - ri.x;
            double ry = r[j].y - ri.y;
            double rz = r[j].z - ri.z;
            double r2 = rx*rx + ry*ry + rz*rz;

            epot_tmp -= (ri.w * r[j].w) * rsqrt(r2);
        }

        double vx = vi.x * vi.x;
        double vy = vi.y * vi.y;
        double vz = vi.z * vi.z;
        double v2 = vx + vy + vz;

        ekin_tmp = 0.5 * ri.w * v2;

        ekin[i] = ekin_tmp;
        epot[i] = epot_tmp;
    }
}

void Hermite4GPU::get_kernel_error(){
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "[Error] : ";
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}

void Hermite4GPU::gpu_timer_start(){
    cudaEventRecord(start);
}

float Hermite4GPU::gpu_timer_stop(std::string f){
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float msec = 0;
    cudaEventElapsedTime(&msec, start, stop);
    #if KERNEL_TIME
    if (f != "")
        std::cout << "[Time] " << f << " : " << msec << " msec" << std::endl;
    #endif
    return msec;
}

void Hermite4GPU::force_calculation(Predictor pi, Predictor pj, Forces &fi)
{
    double rx = pj.r[0] - pi.r[0];
    double ry = pj.r[1] - pi.r[1];
    double rz = pj.r[2] - pi.r[2];

    double vx = pj.v[0] - pi.v[0];
    double vy = pj.v[1] - pi.v[1];
    double vz = pj.v[2] - pi.v[2];

    double r2     = rx*rx + ry*ry + rz*rz + ns->e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * pj.m;
    double mr5inv = r5inv * pj.m;

    double rv = rx*vx + ry*vy + rz*vz;

    fi.a[0] += (rx * mr3inv);
    fi.a[1] += (ry * mr3inv);
    fi.a[2] += (rz * mr3inv);

    fi.a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    fi.a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    fi.a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void Hermite4GPU::integration()
{
    ns->gtime.integration_ini = omp_get_wtime();

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = ns->snapshot_time;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;


    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    init_acc_jrk();
    init_dt(ATIME, ETA_S, ITIME);

    ns->en.ini = get_energy_gpu();   // Initial calculation of the energy of the system
    ns->en.tmp = ns->en.ini;

    //ns->hmr_time = nu->get_half_mass_relaxation_time();
    //ns->cr_time  = nu->get_crossing_time();

    logger->print_info();
    logger->write_info();
    logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, ns->en.ini);

    int snap_number = ns->snapshot_number;
    logger->write_snapshot(snap_number, ITIME);
    snap_number++;

    if (ns->ops.print_all)
    {
        logger->print_all(ITIME);
    }
    if (ns->ops.print_lagrange)
    {
        nu->lagrange_radii();
        logger->print_lagrange_radii(ITIME, nu->layers_radii);
    }

    while (ITIME < ns->integration_time)
    {
        ITIME = ATIME;

        nact = find_particles_to_move(ITIME);

        save_old_acc_jrk(nact);

        predicted_pos_vel(ITIME);

        update_acc_jrk(nact);

        correction_pos_vel(ITIME, nact);

        // Update the amount of interactions counter
        interactions += nact * ns->n;

        // Find the next integration time
        next_integration_time(ATIME);


        if(nact == ns->n)
        {
            //assert(nact == ns->n);
            logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, get_energy_gpu());
            if (ns->ops.print_all)
            {
                logger->print_all(ITIME);
            }
            if (ns->ops.print_lagrange)
            {
                nu->lagrange_radii();
                logger->print_lagrange_radii(ITIME, nu->layers_radii);
            }
            logger->write_snapshot(snap_number, ITIME);
            snap_number++;
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        ns->iterations++;

    }
    ns->gtime.integration_end =  omp_get_wtime() - ns->gtime.integration_ini;
    logger->write_snapshot(snap_number, ITIME);
    //logger->add_info(std::string("SnapshotNumber:"), std::to_string(snap_number));
    logger->add_info(std::string("SnapshotNumber:"), std::string(SSTR(snap_number)));
}
