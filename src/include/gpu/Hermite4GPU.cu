#include "Hermite4GPU.cuh"

void Hermite4GPU::set_pointers(Predictor *dp, Predictor *di, Predictor *hi,
                               Forces *dfout, Forces *dfouttmp, Forces *hfouttmp,
                               Forces *df, int *dmove)
{
    d_p = dp;
    d_i = di;
    h_i = hi;
    d_fout = dfout;
    d_fout_tmp = dfouttmp;
    h_fout_tmp = hfouttmp;
    d_f = df;
    d_move = dmove;
}


void Hermite4GPU::init_acc_jrk(Predictor *p, Forces *f)
{

    CUDA_SAFE_CALL(cudaMemcpy(d_p,
                                p,
                                n * sizeof(Predictor),
                                cudaMemcpyHostToDevice));

    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (d_p, d_f, n, e2);
    //get_kernel_error();

    CUDA_SAFE_CALL(cudaMemcpy(f,
                            d_f,
                            n * sizeof(Forces),
                            cudaMemcpyDeviceToHost));
}

void Hermite4GPU::update_acc_jrk(int nact, int *move, Predictor *p, Forces *f, Gtime &gtime)
{
    gtime.update_ini = omp_get_wtime();

    // Copying to the device the predicted r and v
    CUDA_SAFE_CALL(cudaMemcpy(d_p, p, sizeof(Predictor) * n,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_move, move, sizeof(int) * n,cudaMemcpyHostToDevice));

    // Fill the h_i Predictor array with the particles that we need
    // to move in this iteration
    for (int i = 0; i < nact; i++) {
        int id = move[i];
        h_i[i] = p[id];
    }

    // Copy to the GPU (d_i) the preddictor host array (h_i)
    CUDA_SAFE_CALL(cudaMemcpy(d_i, h_i, sizeof(Predictor) * nact, cudaMemcpyHostToDevice));


    // Blocks, threads and shared memory configuration
    int nact_blocks = 1 + (nact-1)/BSIZE;
    dim3 nblocks(nact_blocks,NJBLOCK, 1);
    dim3 nthreads(BSIZE, 1, 1);

    //std::cout << "nact: " << nact << " , BSIZE: " << BSIZE << std::endl;
    //assert(nact >= BSIZE);
    // Kernel to update the forces for the particles in d_i
    gtime.grav_ini = omp_get_wtime();
    k_update <<< nblocks, nthreads, smem >>> (d_i, d_p, d_fout,d_move,n, nact,e2);
    cudaThreadSynchronize();
    gtime.grav_end += omp_get_wtime() - gtime.grav_ini;
    //get_kernel_error();

    // Blocks, threads and shared memory configuration for the reduction.
    dim3 rgrid   (nact,   1, 1);
    dim3 rthreads(NJBLOCK, 1, 1);

    // Kernel to reduce que temp array with the forces
    gtime.reduce_ini = omp_get_wtime();
    reduce <<< rgrid, rthreads, smem_reduce >>>(d_fout, d_fout_tmp, nact);
    gtime.reduce_end += omp_get_wtime() - gtime.reduce_ini;
    //get_kernel_error();

    // Copy from the GPU the new forces for the d_i particles.
    CUDA_SAFE_CALL(cudaMemcpy(h_fout_tmp, d_fout_tmp, sizeof(Forces) * nact, cudaMemcpyDeviceToHost));

    // Update forces in the host
    for (int i = 0; i < nact; i++) {
        int id = move[i];
        f[id] = h_fout_tmp[i];
    }

    gtime.update_end += (omp_get_wtime() - gtime.update_ini);

}



/*
 * @fn k_init_acc_jrk
 *
 * @desc GPU Kernel which calculates the initial acceleration and jerk
 * of all the particles of the system.
 *
 */
__global__ void k_init_acc_jrk (Predictor *p,
                                Forces *d_f,
                                int n,
                                double e2)
{

    extern __shared__ Predictor sh[];

    Forces ff;
    ff.a[0] = 0.0;
    ff.a[1] = 0.0;
    ff.a[2] = 0.0;
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

        d_f[id] = ff;
    }
}

__device__ void k_force_calculation(Predictor i_p,
                                     Predictor j_p,
                                     Forces &d_f,
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

    d_f.a[0] += (rx * mr3inv);
    d_f.a[1] += (ry * mr3inv);
    d_f.a[2] += (rz * mr3inv);

    d_f.a1[0] += (vx * mr3inv - (3 * rv) * rx * mr5inv);
    d_f.a1[1] += (vy * mr3inv - (3 * rv) * ry * mr5inv);
    d_f.a1[2] += (vz * mr3inv - (3 * rv) * rz * mr5inv);
}

/*
 * @fn k_update()
 *
 * @brief Gravitational interaction kernel.
 */
__global__ void k_update(Predictor *d_i,
                         Predictor *d_j,
                         Forces *d_fout,
                         int *move,
                         int n,
                         int total,
                         double e2)
{
    int ibid = blockIdx.x;
    int jbid = blockIdx.y;
    int tid = threadIdx.x;
    int iaddr = tid + blockDim.x * ibid;
    int jstart = (n * (jbid  )) / NJBLOCK;
    int jend   = (n * (jbid+1)) / NJBLOCK;

    Predictor ip = d_i[iaddr];
    Forces fo;
    fo.a[0] = 0.0;
    fo.a[1] = 0.0;
    fo.a[2] = 0.0;
    fo.a1[0] = 0.0;
    fo.a1[1] = 0.0;
    fo.a1[2] = 0.0;

    //if (iaddr < total || threadIdx.y + jbid * blockDim.y < n)
    //{
        for(int j=jstart; j<jend; j+=BSIZE){
            __shared__ Predictor jpshare[BSIZE];
            __syncthreads();
            Predictor *src = (Predictor *)&d_j[j];
            Predictor *dst = (Predictor *)jpshare;
            dst[      tid] = src[      tid];
            dst[BSIZE+tid] = src[BSIZE+tid];
            __syncthreads();

            // If the total amount of particles is not a multiple of BSIZE
            if(jend-j < BSIZE){
                #pragma unroll 4
                for(int jj=0; jj<jend-j; jj++){
                    Predictor jp = jpshare[jj];
                    k_force_calculation(ip, jp, fo, e2);
                }
            }
            else{
                #pragma unroll 4
                for(int jj=0; jj<BSIZE; jj++){
                    Predictor jp = jpshare[jj];
                    k_force_calculation(ip, jp, fo, e2);
                }
            }
        }

        d_fout[iaddr*NJBLOCK + jbid] = fo;
    //}
}

/*
 * @fn reduce()
 *
 * @brief Forces reduction kernel
 */
__global__ void reduce(Forces *d_in,
                       Forces *d_out,
                       unsigned int total)
{
    extern __shared__ Forces sdata[];

    const int xid   = threadIdx.x;
    const int bid   = blockIdx.x;
    const int iaddr = xid + blockDim.x * bid;

    sdata[xid] = d_in[iaddr];
    __syncthreads();

    if(xid < 8) sdata[xid] += sdata[xid + 8];
    if(xid < 4) sdata[xid] += sdata[xid + 4];
    if(xid < 2) sdata[xid] += sdata[xid + 2];
    if(xid < 1) sdata[xid] += sdata[xid + 1];

    if(xid == 0){
        d_out[bid] = sdata[0];
    }
}


/** Not implemented using GPU */
void Hermite4GPU::predicted_pos_vel(double ITIME, Predictor *p, double4 *r, double4 *v, Forces *f, double *t, Gtime &gtime)
{

    gtime.prediction_ini = omp_get_wtime();
    for (int i = 0; i < n; i++)
    {
        double dt  = ITIME - t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        p[i].r[0] = (dt3/6 * f[i].a1[0]) + (dt2/2 * f[i].a[0]) + (dt * v[i].x) + r[i].x;
        p[i].r[1] = (dt3/6 * f[i].a1[1]) + (dt2/2 * f[i].a[1]) + (dt * v[i].y) + r[i].y;
        p[i].r[2] = (dt3/6 * f[i].a1[2]) + (dt2/2 * f[i].a[2]) + (dt * v[i].z) + r[i].z;

        p[i].v[0] = (dt2/2 * f[i].a1[0]) + (dt * f[i].a[0]) + v[i].x;
        p[i].v[1] = (dt2/2 * f[i].a1[1]) + (dt * f[i].a[1]) + v[i].y;
        p[i].v[2] = (dt2/2 * f[i].a1[2]) + (dt * f[i].a[2]) + v[i].z;

        p[i].m = r[i].w;

    }
    gtime.prediction_end += omp_get_wtime() - gtime.prediction_ini;
}

/** Not implemented using GPU */
void Hermite4GPU::correction_pos_vel(double ITIME, int nact, int *move, double4 *r, double4 *v, Forces *f, double *t, double *dt, Predictor *p, Forces *old, double4 *a3, double4 *a2, Gtime &gtime)
{
    gtime.correction_ini = omp_get_wtime();
    for (int k = 0; k < nact; k++)
    {
        int i = move[k];

        double dt1 = dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        a2[i].x = (-6 * (old[i].a[0] - f[i].a[0] ) - dt1 * (4 * old[i].a1[0] + 2 * f[i].a1[0]) ) / dt2;
        a2[i].y = (-6 * (old[i].a[1] - f[i].a[1] ) - dt1 * (4 * old[i].a1[1] + 2 * f[i].a1[1]) ) / dt2;
        a2[i].z = (-6 * (old[i].a[2] - f[i].a[2] ) - dt1 * (4 * old[i].a1[2] + 2 * f[i].a1[2]) ) / dt2;

        // Acceleration 3rd derivate
        a3[i].x = (12 * (old[i].a[0] - f[i].a[0] ) + 6 * dt1 * (old[i].a1[0] + f[i].a1[0]) ) / dt3;
        a3[i].y = (12 * (old[i].a[1] - f[i].a[1] ) + 6 * dt1 * (old[i].a1[1] + f[i].a1[1]) ) / dt3;
        a3[i].z = (12 * (old[i].a[2] - f[i].a[2] ) + 6 * dt1 * (old[i].a1[2] + f[i].a1[2]) ) / dt3;

        // Correcting position
        r[i].x = p[i].r[0] + (dt4/24)*a2[i].x + (dt5/120)*a3[i].x;
        r[i].y = p[i].r[1] + (dt4/24)*a2[i].y + (dt5/120)*a3[i].y;
        r[i].z = p[i].r[2] + (dt4/24)*a2[i].z + (dt5/120)*a3[i].z;

        // Correcting velocity
        v[i].x = p[i].v[0] + (dt3/6)*a2[i].x + (dt4/24)*a3[i].x;
        v[i].y = p[i].v[1] + (dt3/6)*a2[i].y + (dt4/24)*a3[i].y;
        v[i].z = p[i].v[2] + (dt3/6)*a2[i].z + (dt4/24)*a3[i].z;


        t[i] = ITIME;
        double normal_dt  = get_timestep_normal(i, a2, a3, dt, f, eta);
        normal_dt = normalize_dt(normal_dt, dt[i], t[i], i);
        dt[i] = normal_dt;

    }
    gtime.correction_end += omp_get_wtime() - gtime.correction_ini;
}

