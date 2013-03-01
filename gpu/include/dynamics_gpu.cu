#include "dynamics_gpu.cuh"

inline __host__ __device__ Forces operator+(volatile Forces &a, volatile Forces &b)
{
    Forces tmp;
    tmp.a[0] = a.a[0] + b.a[0];
    tmp.a[1] = a.a[1] + b.a[1];
    tmp.a[2] = a.a[2] + b.a[2];

    tmp.a1[0] = a.a1[0] + b.a1[0];
    tmp.a1[1] = a.a1[1] + b.a1[1];
    tmp.a1[2] = a.a1[2] + b.a1[2];

    return tmp;
}

inline __host__ __device__ void operator+=(volatile Forces &a, Forces &b)
{
    a.a[0] += b.a[0];
    a.a[1] += b.a[1];
    a.a[2] += b.a[2];

    a.a1[0] += b.a1[0];
    a.a1[1] += b.a1[1];
    a.a1[2] += b.a1[2];
}

inline __host__ __device__ void operator+=(volatile Forces &a, volatile Forces &b)
{
    a.a[0] += b.a[0];
    a.a[1] += b.a[1];
    a.a[2] += b.a[2];

    a.a1[0] += b.a1[0];
    a.a1[1] += b.a1[1];
    a.a1[2] += b.a1[2];
}


__host__ void gpu_init_acc_jrk()
{
    int smem = BSIZE * 2* sizeof(double4);
    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (d_r, d_v, d_f, d_m, n);
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_init_acc_jrk: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}

__host__ double gpu_energy()
{

    CUDA_SAFE_CALL(cudaMemcpy(d_r,  h_r,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v,  h_v,  d4_size,cudaMemcpyHostToDevice));
    k_energy <<< nblocks, nthreads >>> (d_r, d_v, d_ekin, d_epot, d_m, n);
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_energy: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif

    CUDA_SAFE_CALL(cudaMemcpy(h_ekin, d_ekin, d1_size,cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_epot, d_epot, d1_size,cudaMemcpyDeviceToHost));

    // Reduction on CPU
    ekin = 0.0;
    epot = 0.0;

    for (int i = 0; i < n; i++) {
        ekin += h_ekin[i];
        epot += h_epot[i];
    }
    return ekin + epot;
}

__host__ void gpu_update_acc_jrk_simple(int total) {

    int smem = BSIZE * sizeof(Predictor);
    //Predictor tmp[total];

    //dim3 nthreads(BSIZE, 1, 1);
    //dim3 nblocks(1 + (total-1)/BSIZE,NJBLOCK, 1);
    k_update_acc_jrk_simple <<< nblocks, nthreads, smem >>> (d_p, d_f,
                                                             d_m, d_move, n, total);
    //cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_update_acc_jrk_simple: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}


/*
#define NXREDUCE 16 // must be >NJBLOCK
#define NYREDUCE  8
 *
__global__ void force_reduce_kernel(
        const int ni,
        const Force fpart[][NJBLOCK],
        __out Force ftot []){
    const int xid = threadIdx.x;
    const int yid = threadIdx.y;
    const int bid = blockIdx.x;
    const int iaddr = yid + blockDim.y * bid;

    __shared__ Force fshare[NYREDUCE][NXREDUCE];
    if(xid < NJBLOCK){
        fshare[yid][xid] = fpart[iaddr][xid];
    }else{
        fshare[yid][xid].clear();
    }
    Force *fs = fshare[yid];
#if NXREDUCE==32
    if(xid < 16) fs[xid] += fs[xid + 16];
#endif
    if(xid < 8) fs[xid] += fs[xid + 8];
    if(xid < 4) fs[xid] += fs[xid + 4];
    if(xid < 2) fs[xid] += fs[xid + 2];
    if(xid < 1) fs[xid] += fs[xid + 1];

    if(iaddr < ni){
        ftot[iaddr] = fs[0];
    }
}
*/

#define NXREDUCE 16
#define NYREDUCE  8

__global__ void reduce(Forces *d_in, Forces *d_out, unsigned int n)
{
    __shared__ Forces sdata[NYREDUCE * NXREDUCE];

    const int xid = threadIdx.x;
    const int yid = threadIdx.y;
    const int bid = blockIdx.x;
    const int iaddr = yid + blockDim.y * bid;

    if(xid < NJBLOCK){
        //fshare[yid][xid] = fpart[iaddr][xid];
        sdata[yid*NYREDUCE + xid] = d_in[iaddr * NYREDUCE + xid];
    }else{
        //fshare[yid][xid].clear();
        sdata[yid*NYREDUCE +xid].a[0] = 0.0;
        sdata[yid*NYREDUCE +xid].a[1] = 0.0;
        sdata[yid*NYREDUCE +xid].a[2] = 0.0;
        sdata[yid*NYREDUCE +xid].a1[0] = 0.0;
        sdata[yid*NYREDUCE +xid].a1[1] = 0.0;
        sdata[yid*NYREDUCE +xid].a1[2] = 0.0;
    }
    //Forces *fs = fshare[yid];
    Forces *fs = sdata;
//#if NXREDUCE==32
//    if(xid < 16) fs[xid] += fs[xid + 16];
//#endif
    if(xid < 8) fs[yid * NYREDUCE + xid] += fs[yid * NYREDUCE + xid + 8];
    if(xid < 4) fs[yid * NYREDUCE + xid] += fs[yid * NYREDUCE + xid + 4];
    if(xid < 2) fs[yid * NYREDUCE + xid] += fs[yid * NYREDUCE + xid + 2];
    if(xid < 1) fs[yid * NYREDUCE + xid] += fs[yid * NYREDUCE + xid + 1];

    if(iaddr < n){
        d_out[iaddr] = fs[0];
    }

}


__host__ void gpu_update(int total) {

    // Fill the h_i Predictor array with the particles that we need
    // to move in this iteration
    for (int i = 0; i < total; i++) {
        int id = h_move[i];
        h_i[i] = h_p[id];
    }

    // Copy to the GPU (d_i) the preddictor host array (h_i)
    CUDA_SAFE_CALL(cudaMemcpy(d_i, h_i, sizeof(Predictor) * total, cudaMemcpyHostToDevice));


    // Blocks and Threads configuration
    dim3 nblocks2(1 + (total-1)/BSIZE,NJBLOCK, 1);
    dim3 nthreads2(BSIZE, 1, 1);
    int smem = BSIZE * sizeof(Predictor);

    // Kernel call
    k_update <<< nblocks2, nthreads2, smem >>> (d_i, d_p, d_fout,d_m, n, total);
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_update: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif


    // cpu
    //float t1 = gettime_ms;
    CUDA_SAFE_CALL(cudaMemcpy(h_fout, d_fout, sizeof(Forces) * total * NJBLOCK, cudaMemcpyDeviceToHost));
    // Reduction
    Forces tmp_f, new_f;
    int id;
    for (int i = 0; i < total ; i++) {
        id = h_move[i];
        memset(&tmp_f, 0, sizeof(Forces));
        for (int j = 0; j < NJBLOCK; j++) {
            new_f = h_fout[i * NJBLOCK + j];
            tmp_f  += new_f;
            //tmp_f.a[0]  += new_f.a[0];
            //tmp_f.a[1]  += new_f.a[1];
            //tmp_f.a[2]  += new_f.a[2];
            //tmp_f.a1[0] += new_f.a1[0];
            //tmp_f.a1[1] += new_f.a1[1];
            //tmp_f.a1[2] += new_f.a1[2];
        }
        h_f[id] = tmp_f;
     //   printf("%f ", h_f[i].a[0]);
    }
    //printf("\n");
    //float t2 = gettime_ms;
    //printf("Reduction CPU: %f\n", t2-t1);
    // end cpu


    //// begin gpu
    //t1 = gettime_ms;
    ////size_t nb = n/NJBLOCK;
    ////reduce <<< nb, NJBLOCK, smem >>>(d_fout, d_fout_tmp, n);
    //const int ni8 = 1 + (total-1) / NYREDUCE;
    //dim3 rgrid   (ni8, 1, 1);
    //dim3 rthreads(NXREDUCE, NYREDUCE, 1);
    //reduce <<< rgrid, rthreads >>>(d_fout, d_fout_tmp, total);
    //CUDA_SAFE_CALL(cudaMemcpy(h_fout_tmp, d_fout_tmp, sizeof(Forces) * total * NJBLOCK, cudaMemcpyDeviceToHost));
    //for (int i = 0; i < total; i++) {
    //    printf("%f ", h_fout_tmp[i].a[0]);
    //}
    //printf("\n");
    //t2 = gettime_ms;
    //printf("Reduction GPU: %f\n", t2-t1);


    //// end gpu

    //getchar();

}
