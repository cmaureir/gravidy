#include "dynamics_gpu.cuh"


/*
 * @fn __host__ void gpu_init_acc_jrk()
 *
 * @brief
 *  Initial calculation of the acceleration and jerk on the GPU.
 */
__host__ void gpu_init_acc_jrk()
{
    int smem = BSIZE * 2* sizeof(double4);

    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (d_r, d_v, d_f, d_m, n,e2);
    cudaThreadSynchronize();
    get_kernel_error();
}

/*
 * @fn __host__ double gpu_energy()
 *
 * @brief
 *  Energy calculation on the GPU.
 */
__host__ double gpu_energy()
{

    CUDA_SAFE_CALL(cudaMemcpy(d_r,  h_r,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v,  h_v,  d4_size,cudaMemcpyHostToDevice));

    gpu_timer_start();
    k_energy <<< nblocks, nthreads >>> (d_r, d_v, d_ekin, d_epot, d_m, n);
    cudaThreadSynchronize();
    float msec = gpu_timer_stop("k_energy");
    get_kernel_error();

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

/*
 * @fn __host__ void gpu_update()
 *
 * @brief
 *  Gravitational interactions calculation for the nact particles
 *  using j-parallelization, and then doing a reduction.
 */
__host__ void gpu_update(int total)
{

    gtime.update_ini = omp_get_wtime();

    // Copying to the device the predicted r and v
    CUDA_SAFE_CALL(cudaMemcpy(d_p, h_p, sizeof(Predictor) * n,cudaMemcpyHostToDevice));

    // Fill the h_i Predictor array with the particles that we need
    // to move in this iteration
    for (int i = 0; i < total; i++) {
        int id = h_move[i];
        h_i[i] = h_p[id];
    }

    // Copy to the GPU (d_i) the preddictor host array (h_i)
    CUDA_SAFE_CALL(cudaMemcpy(d_i, h_i, sizeof(Predictor) * total, cudaMemcpyHostToDevice));


    // Blocks, threads and shared memory configuration
    dim3 nblocks(1 + (total-1)/BSIZE,NJBLOCK, 1);
    dim3 nthreads(BSIZE, 1, 1);
    size_t smem = BSIZE * sizeof(Predictor);

    // Kernel to update the forces for the particles in d_i
    gtime.grav_ini = omp_get_wtime();
    k_update <<< nblocks, nthreads, smem >>> (d_i, d_p, d_fout,d_m, n, total,e2);
    cudaThreadSynchronize();
    gtime.grav_end += omp_get_wtime() - gtime.grav_ini;
    get_kernel_error();

    // Blocks, threads and shared memory configuration for the reduction.
    dim3 rgrid   (total,   1, 1);
    dim3 rthreads(NJBLOCK, 1, 1);
    size_t smem2 = sizeof(Forces) * NJBLOCK + 1;

    // Kernel to reduce que temp array with the forces
    gtime.reduce_ini = omp_get_wtime();
    reduce <<< rgrid, rthreads, smem2 >>>(d_fout, d_fout_tmp, total);
    cudaThreadSynchronize();
    gtime.reduce_end += omp_get_wtime() - gtime.grav_ini;
    get_kernel_error();

    // Copy from the GPU the new forces for the d_i particles.
    CUDA_SAFE_CALL(cudaMemcpy(h_fout_tmp, d_fout_tmp, sizeof(Forces) * total, cudaMemcpyDeviceToHost));

    // Update forces in the host
    for (int i = 0; i < total; i++) {
        int id = h_move[i];
        h_f[id] = h_fout_tmp[i];
    }

    gtime.update_end += omp_get_wtime() - gtime.update_ini;

}
