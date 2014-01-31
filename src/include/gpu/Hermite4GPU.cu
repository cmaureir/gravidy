#include "Hermite4GPU.cuh"

Hermite4GPU::~Hermite4GPU()
{
    free_arrays_device();
}

void Hermite4GPU::alloc_arrays_device()
{
    int d4_size = ns->n * sizeof(double4);
    int d1_size = ns->n * sizeof(double);
    int i1_size = ns->n * sizeof(int);
    int ff_size = ns->n * sizeof(Forces);
    int pp_size = ns->n * sizeof(Predictor);

    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_r,        d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_v,        d4_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_f,        ff_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_p,        pp_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_ekin,     d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_epot,     d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_t,        d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_dt,       d1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_move,     i1_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_i,        pp_size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_fout,     ff_size * NJBLOCK));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ns->d_fout_tmp, ff_size * NJBLOCK));

    CUDA_SAFE_CALL(cudaMemset(ns->d_r,         0, d4_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_v,         0, d4_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_f,         0, ff_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_p,         0, pp_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_ekin,      0, d1_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_epot,      0, d1_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_t,         0, d1_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_dt,        0, d1_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_move,      0, i1_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_i,         0, pp_size));
    CUDA_SAFE_CALL(cudaMemset(ns->d_fout,      0, ff_size * NJBLOCK));
    CUDA_SAFE_CALL(cudaMemset(ns->d_fout_tmp,  0, ff_size * NJBLOCK));

    // Extra CPU array
    ns->h_fout_tmp= new Forces[ff_size*NJBLOCK];

}

void Hermite4GPU::free_arrays_device()
{
    CUDA_SAFE_CALL(cudaFree(ns->d_r));
    CUDA_SAFE_CALL(cudaFree(ns->d_v));
    CUDA_SAFE_CALL(cudaFree(ns->d_f));
    CUDA_SAFE_CALL(cudaFree(ns->d_p));
    CUDA_SAFE_CALL(cudaFree(ns->d_ekin));
    CUDA_SAFE_CALL(cudaFree(ns->d_epot));
    CUDA_SAFE_CALL(cudaFree(ns->d_t));
    CUDA_SAFE_CALL(cudaFree(ns->d_dt));
    CUDA_SAFE_CALL(cudaFree(ns->d_move));
    CUDA_SAFE_CALL(cudaFree(ns->d_i));
    CUDA_SAFE_CALL(cudaFree(ns->d_fout));
    CUDA_SAFE_CALL(cudaFree(ns->d_fout_tmp));

    delete ns->h_fout_tmp;
}

/** Not implemented using GPU */
void Hermite4GPU::predicted_pos_vel(double ITIME)
{

    ns->gtime.prediction_ini = omp_get_wtime();
    for (int i = 0; i < ns->n; i++)
    {
        double dt  = ITIME - ns->h_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        ns->h_p[i].r[0] = (dt3/6 * ns->h_f[i].a1[0]) + (dt2/2 * ns->h_f[i].a[0]) + (dt * ns->h_v[i].x) + ns->h_r[i].x;
        ns->h_p[i].r[1] = (dt3/6 * ns->h_f[i].a1[1]) + (dt2/2 * ns->h_f[i].a[1]) + (dt * ns->h_v[i].y) + ns->h_r[i].y;
        ns->h_p[i].r[2] = (dt3/6 * ns->h_f[i].a1[2]) + (dt2/2 * ns->h_f[i].a[2]) + (dt * ns->h_v[i].z) + ns->h_r[i].z;

        ns->h_p[i].v[0] = (dt2/2 * ns->h_f[i].a1[0]) + (dt * ns->h_f[i].a[0]) + ns->h_v[i].x;
        ns->h_p[i].v[1] = (dt2/2 * ns->h_f[i].a1[1]) + (dt * ns->h_f[i].a[1]) + ns->h_v[i].y;
        ns->h_p[i].v[2] = (dt2/2 * ns->h_f[i].a1[2]) + (dt * ns->h_f[i].a[2]) + ns->h_v[i].z;

        ns->h_p[i].m = ns->h_r[i].w;

    }
    ns->gtime.prediction_end += omp_get_wtime() - ns->gtime.prediction_ini;
}

/** Not implemented using GPU */
void Hermite4GPU::correction_pos_vel(double ITIME, int nact)
{
    ns->gtime.correction_ini = omp_get_wtime();
    for (int k = 0; k < nact; k++)
    {
        int i = ns->h_move[k];

        double dt1 = ns->h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        ns->h_a2[i].x = (-6 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) - dt1 * (4 * ns->h_old[i].a1[0] + 2 * ns->h_f[i].a1[0]) ) / dt2;
        ns->h_a2[i].y = (-6 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) - dt1 * (4 * ns->h_old[i].a1[1] + 2 * ns->h_f[i].a1[1]) ) / dt2;
        ns->h_a2[i].z = (-6 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) - dt1 * (4 * ns->h_old[i].a1[2] + 2 * ns->h_f[i].a1[2]) ) / dt2;

        // Acceleration 3rd derivate
        ns->h_a3[i].x = (12 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) + 6 * dt1 * (ns->h_old[i].a1[0] + ns->h_f[i].a1[0]) ) / dt3;
        ns->h_a3[i].y = (12 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) + 6 * dt1 * (ns->h_old[i].a1[1] + ns->h_f[i].a1[1]) ) / dt3;
        ns->h_a3[i].z = (12 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) + 6 * dt1 * (ns->h_old[i].a1[2] + ns->h_f[i].a1[2]) ) / dt3;

        // Correcting position
        ns->h_r[i].x = ns->h_p[i].r[0] + (dt4/24)*ns->h_a2[i].x + (dt5/120)*ns->h_a3[i].x;
        ns->h_r[i].y = ns->h_p[i].r[1] + (dt4/24)*ns->h_a2[i].y + (dt5/120)*ns->h_a3[i].y;
        ns->h_r[i].z = ns->h_p[i].r[2] + (dt4/24)*ns->h_a2[i].z + (dt5/120)*ns->h_a3[i].z;

        // Correcting velocity
        ns->h_v[i].x = ns->h_p[i].v[0] + (dt3/6)*ns->h_a2[i].x + (dt4/24)*ns->h_a3[i].x;
        ns->h_v[i].y = ns->h_p[i].v[1] + (dt3/6)*ns->h_a2[i].y + (dt4/24)*ns->h_a3[i].y;
        ns->h_v[i].z = ns->h_p[i].v[2] + (dt3/6)*ns->h_a2[i].z + (dt4/24)*ns->h_a3[i].z;


        ns->h_t[i] = ITIME;
        double normal_dt  = nu->get_timestep_normal(i);
        normal_dt = nu->normalize_dt(normal_dt, ns->h_dt[i], ns->h_t[i], i);
        ns->h_dt[i] = normal_dt;

    }
    ns->gtime.correction_end += omp_get_wtime() - ns->gtime.correction_ini;
}

void Hermite4GPU::init_acc_jrk()
{

    CUDA_SAFE_CALL(cudaMemcpy(ns->d_p,
                              ns->h_p,
                              ns->n * sizeof(Predictor),
                              cudaMemcpyHostToDevice));

    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (ns->d_p,
                                                    ns->d_f,
                                                    ns->n,
                                                    ns->e2);
    //get_kernel_error();

    CUDA_SAFE_CALL(cudaMemcpy(ns->h_f,
                              ns->d_f,
                              ns->n * sizeof(Forces),
                              cudaMemcpyDeviceToHost));
}

void Hermite4GPU::update_acc_jrk(int nact)
{
    ns->gtime.update_ini = omp_get_wtime();

    // Copying to the device the predicted r and v
    CUDA_SAFE_CALL(cudaMemcpy(ns->d_p,
                              ns->h_p,
                              ns->n * sizeof(Predictor),
                              cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(ns->d_move,
                              ns->h_move,
                              ns->n * sizeof(int),
                              cudaMemcpyHostToDevice));

    // Fill the h_i Predictor array with the particles that we need
    // to move in this iteration
    for (int i = 0; i < nact; i++)
    {
        int id = ns->h_move[i];
        ns->h_i[i] = ns->h_p[id];
    }

    // Copy to the GPU (d_i) the preddictor host array (h_i)
    CUDA_SAFE_CALL(cudaMemcpy(ns->d_i,
                              ns->h_i,
                              nact * sizeof(Predictor),
                              cudaMemcpyHostToDevice));


    // Blocks, threads and shared memory configuration
    int  nact_blocks = 1 + (nact-1)/BSIZE;
    dim3 nblocks(nact_blocks,NJBLOCK, 1);
    dim3 nthreads(BSIZE, 1, 1);

    // Kernel to update the forces for the particles in d_i
    ns->gtime.grav_ini = omp_get_wtime();
    k_update <<< nblocks, nthreads, smem >>> (ns->d_i,
                                              ns->d_p,
                                              ns->d_fout,
                                              ns->d_move,
                                              ns->n,
                                              nact,
                                              ns->e2);
    cudaThreadSynchronize();
    ns->gtime.grav_end += omp_get_wtime() - ns->gtime.grav_ini;
    //get_kernel_error();

    // Blocks, threads and shared memory configuration for the reduction.
    dim3 rgrid   (nact,   1, 1);
    dim3 rthreads(NJBLOCK, 1, 1);

    // Kernel to reduce que temp array with the forces
    ns->gtime.reduce_ini = omp_get_wtime();
    reduce <<< rgrid, rthreads, smem_reduce >>>(ns->d_fout,
                                                ns->d_fout_tmp);
    ns->gtime.reduce_end += omp_get_wtime() - ns->gtime.reduce_ini;
    //get_kernel_error();

    // Copy from the GPU the new forces for the d_i particles.
    CUDA_SAFE_CALL(cudaMemcpy(ns->h_fout_tmp,
                              ns->d_fout_tmp,
                              nact * sizeof(Forces),
                              cudaMemcpyDeviceToHost));

    // Update forces in the host
    for (int i = 0; i < nact; i++)
    {
        int id = ns->h_move[i];
        ns->h_f[id] = ns->h_fout_tmp[i];
    }

    ns->gtime.update_end += (omp_get_wtime() - ns->gtime.update_ini);
}

double Hermite4GPU::get_energy_gpu()
{
    CUDA_SAFE_CALL(cudaMemcpy(ns->d_r, ns->h_r,  sizeof(double4) * ns->n,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(ns->d_v, ns->h_v,  sizeof(double4) * ns->n,cudaMemcpyHostToDevice));

    //gpu_timer_start();
    int nthreads = BSIZE;
    int nblocks = std::ceil(ns->n/(float)nthreads);
    k_energy <<< nblocks, nthreads >>> (ns->d_r,
                                        ns->d_v,
                                        ns->d_ekin,
                                        ns->d_epot,
                                        ns->n,
                                        ns->e2);
    cudaThreadSynchronize();
    //float msec = gpu_timer_stop("k_energy");
    //get_kernel_error();

    CUDA_SAFE_CALL(cudaMemcpy(ns->h_ekin, ns->d_ekin, sizeof(double) * ns->n,cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(ns->h_epot, ns->d_epot, sizeof(double) * ns->n,cudaMemcpyDeviceToHost));

    // Reduction on CPU
    ns->en.kinetic = 0.0;
    ns->en.potential = 0.0;

    for (int i = 0; i < ns->n; i++)
    {
        ns->en.kinetic   += ns->h_ekin[i];
        ns->en.potential += ns->h_epot[i];
    }
    return ns->en.kinetic + ns->en.potential;
}

/*
 * @fn k_init_acc_jrk
 *
 * @desc GPU Kernel which calculates the initial acceleration and jerk
 * of all the particles of the system.
 *
 */
__global__ void k_init_acc_jrk (Predictor *p,
                                Forces *f,
                                int n,
                                double e2)
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

    if (id < n)
    {
        Predictor pred = p[id];
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
                         int *move,
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
                         double e2)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j;
    double ekin_tmp = 0.0;

    if (i < n)
    {
        double epot_tmp = 0.0;
        for (j = i+1; j < n; j++)
        {
            double rx = r[j].x - r[i].x;
            double ry = r[j].y - r[i].y;
            double rz = r[j].z - r[i].z;
            double r2 = rx*rx + ry*ry + rz*rz + e2;

            epot_tmp -= (r[i].w * r[j].w) * rsqrt(r2);
        }

        double vx = v[i].x * v[i].x;
        double vy = v[i].y * v[i].y;
        double vz = v[i].z * v[i].z;
        double v2 = vx + vy + vz;

        ekin_tmp = 0.5 * r[i].w * v2;

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

void Hermite4GPU::force_calculation(int i, int j)
{
    double rx = ns->h_p[j].r[0] - ns->h_p[i].r[0];
    double ry = ns->h_p[j].r[1] - ns->h_p[i].r[1];
    double rz = ns->h_p[j].r[2] - ns->h_p[i].r[2];

    double vx = ns->h_p[j].v[0] - ns->h_p[i].v[0];
    double vy = ns->h_p[j].v[1] - ns->h_p[i].v[1];
    double vz = ns->h_p[j].v[2] - ns->h_p[i].v[2];

    double r2     = rx*rx + ry*ry + rz*rz + ns->e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * ns->h_p[j].m;
    double mr5inv = r5inv * ns->h_p[j].m;

    double rv = rx*vx + ry*vy + rz*vz;

    ns->h_f[i].a[0] += (rx * mr3inv);
    ns->h_f[i].a[1] += (ry * mr3inv);
    ns->h_f[i].a[2] += (rz * mr3inv);

    ns->h_f[i].a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    ns->h_f[i].a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    ns->h_f[i].a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void Hermite4GPU::integration()
{
    ns->gtime.integration_ini = omp_get_wtime();

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;


    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    init_acc_jrk();
    init_dt(ATIME);

    ns->en.ini = get_energy_gpu();   // Initial calculation of the energy of the system
    ns->en.tmp = ns->en.ini;

    //ns->hmr_time = nu->get_half_mass_relaxation_time();
    //ns->cr_time  = nu->get_crossing_time();

    logger->print_info();
    logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, ns->en.ini);

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


        if(std::ceil(ITIME) == ITIME)
        {
            assert(nact == ns->n);
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
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        ns->iterations++;
    }
    ns->gtime.integration_end =  omp_get_wtime() - ns->gtime.integration_ini;

}
