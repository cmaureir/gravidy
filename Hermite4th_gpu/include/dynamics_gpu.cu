#include "dynamics_gpu.cuh"

__host__ void gpu_init_acc_jrk()
{
    int smem = BSIZE * 2* sizeof(double4);
    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (d_r, d_v, d_f, d_m, n);
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_init_acc_jrk: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}

__host__ void gpu_get_energy_log(double ITIME, int iterations, int nsteps)
{
    energy_end = gpu_energy();
    double relative_error = abs((energy_end-energy_ini)/energy_ini);
    energy_tmp += relative_error;

    if((int)ITIME == 1)
        printf("ITIME  ITER NSTEPS/ITER    ETA     TIME            ENERGY     REL_ENER     CUM_ENER\n");

    printf("%4d %6d %.6f %f %.6f %.15e %.15e %.15e\n",
            (int)ITIME,
            iterations,
            nsteps/(float)iterations,
            ETA_N,
            (float)clock()/CLOCKS_PER_SEC - ini_time,
            energy_end,
            relative_error,
            energy_tmp/(int)ITIME);
}

__host__ double gpu_energy()
{
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

__host__ void gpu_update(int total) {

    int smem = BSIZE * sizeof(Predictor);
    for (int i = 0; i < total; i++) {
        int id = h_move[i];
        h_i[i] = h_p[id];
    }
    CUDA_SAFE_CALL(cudaMemcpy(d_i, h_i, sizeof(Predictor) * total, cudaMemcpyHostToDevice));
    dim3 nblocks2(1 + (total-1)/BSIZE,NJBLOCK, 1);
    dim3 nthreads2(BSIZE, 1, 1);
    k_update <<< nblocks2, nthreads2, smem >>> (d_i, d_p, d_fout,d_m, n, total);
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_update: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
    CUDA_SAFE_CALL(cudaMemcpy(h_fout, d_fout, sizeof(Forces) * total * NJBLOCK, cudaMemcpyDeviceToHost));
    //CUDA_SAFE_CALL(cudaMemcpy(h_fout, d_fout, sizeof(Forces) * n * NJBLOCK, cudaMemcpyDeviceToHost));

    for (int i = 0; i < total * NJBLOCK; i++) {
        printf("%f ", h_fout[i].a[0]);
        if (i%NJBLOCK == 0) printf("||\n");
    }
    getchar();
    // Reduction
    Forces tmp_f;
    for (int i = 0; i < total; i++) {
        int id = h_move[i];
        tmp_f.a[0] = 0.0;
        tmp_f.a[1] = 0.0;
        tmp_f.a[2] = 0.0;
        tmp_f.a1[0] = 0.0;
        tmp_f.a1[1] = 0.0;
        tmp_f.a1[2] = 0.0;
        for (int j = 0; j < NJBLOCK; j++) {
            tmp_f.a[0]  += h_fout[i * NJBLOCK + j].a[0];
            tmp_f.a[1]  += h_fout[i * NJBLOCK + j].a[1];
            tmp_f.a[2]  += h_fout[i * NJBLOCK + j].a[2];
            tmp_f.a1[0] += h_fout[i * NJBLOCK + j].a1[0];
            tmp_f.a1[1] += h_fout[i * NJBLOCK + j].a1[1];
            tmp_f.a1[2] += h_fout[i * NJBLOCK + j].a1[2];
        }
        h_f[id] = tmp_f;
    printf("Updating %d - %f\t%f\t%f - %f\t%f\t%f\n", id, h_f[id].a[0], h_f[id].a[1], h_f[id].a[2], h_f[id].a1[0], h_f[id].a1[1], h_f[id].a1[2]);
    getchar();
    }

}
