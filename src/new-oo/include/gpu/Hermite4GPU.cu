#include "Hermite4GPU.cuh"

void Hermite4GPU::set_pointers(Predictor *dp, Predictor *di, Predictor *hi,
                               Forces *dfout, Forces *dfouttmp, Forces *hfouttmp)
{
    d_p = dp;
    d_i = di;
    h_i = hi;
    d_fout = dfout;
    d_fout_tmp = dfouttmp;
    h_fout_tmp = hfouttmp;
}


void Hermite4GPU::init_acc_jrk(Predictor *p, Forces *f)
{

    int smem = BSIZE * sizeof(Predictor);
    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (p, f, n, e2);
    cudaThreadSynchronize();
}

void Hermite4GPU::update_acc_jrk(int nact, int *move, Predictor *p, Forces *f, Gtime &gtime)
{
    gtime.update_ini = omp_get_wtime();

    // Copying to the device the predicted r and v
    CUDA_SAFE_CALL(cudaMemcpy(d_p, p, sizeof(Predictor) * n,cudaMemcpyHostToDevice));

    // Fill the h_i Predictor array with the particles that we need
    // to move in this iteration
    for (int i = 0; i < nact; i++) {
        int id = move[i];
        h_i[i] = p[id];
    }

    // Copy to the GPU (d_i) the preddictor host array (h_i)
    CUDA_SAFE_CALL(cudaMemcpy(d_i, h_i, sizeof(Predictor) * nact, cudaMemcpyHostToDevice));


    // Blocks, threads and shared memory configuration
    dim3 nblocks(1 + (nact-1)/BSIZE,NJBLOCK, 1);
    dim3 nthreads(BSIZE, 1, 1);
    size_t smem = BSIZE * sizeof(Predictor);

    // Kernel to update the forces for the particles in d_i
    gtime.grav_ini = omp_get_wtime();
    k_update <<< nblocks, nthreads, smem >>> (d_i, d_p, d_fout,n, nact,e2);
    cudaThreadSynchronize();
    gtime.grav_end += omp_get_wtime() - gtime.grav_ini;
    //get_kernel_error();

    // Blocks, threads and shared memory configuration for the reduction.
    dim3 rgrid   (nact,   1, 1);
    dim3 rthreads(NJBLOCK, 1, 1);
    size_t smem2 = sizeof(Forces) * NJBLOCK + 1;

    // Kernel to reduce que temp array with the forces
    gtime.reduce_ini = omp_get_wtime();
    reduce <<< rgrid, rthreads, smem2 >>>(d_fout, d_fout_tmp, nact);
    cudaThreadSynchronize();
    gtime.reduce_end += omp_get_wtime() - gtime.grav_ini;
    //get_kernel_error();

    // Copy from the GPU the new forces for the d_i particles.
    CUDA_SAFE_CALL(cudaMemcpy(h_fout_tmp, d_fout_tmp, sizeof(Forces) * nact, cudaMemcpyDeviceToHost));

    // Update forces in the host
    for (int i = 0; i < nact; i++) {
        int id = move[i];
        f[id] = h_fout_tmp[i];
    }

    gtime.update_end += omp_get_wtime() - gtime.update_ini;

}


/*
 * @fn k_energy
 *
 * @desc GPU Kernel which calculates the energy of the system.
 *
 */
__global__ void k_energy(double4 *r,
                         double4 *v,
                         double *ekin,
                         double *epot,
                         float *m,
                         int n)
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
    //Predictor *sp = (Predictor*)sh;

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
            //#pragma unroll 4
            for(int jj=0; jj<jend-j; jj++){
                Predictor jp = jpshare[jj];
                k_force_calculation(ip, jp, fo, e2);
            }
        }
        else{
            //#pragma unroll 4
            for(int jj=0; jj<BSIZE; jj++){
                Predictor jp = jpshare[jj];
                k_force_calculation(ip, jp, fo, e2);
            }
        }
    }
    d_fout[iaddr*NJBLOCK + jbid] = fo;

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
