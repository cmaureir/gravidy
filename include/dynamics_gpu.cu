#include "dynamics_gpu.cuh"



__host__ void gpu_update_2d(int total)
{
    // i-threads = number of particles to move (total)
    // j-threads = particles per block to calculate force on i (~ 32)
    ////////ENVIAR SOLO PARTICULAS QUE SE MUEVEN
    cudaMemcpy(d_move, h_move, n * sizeof(int), cudaMemcpyHostToDevice);

    printf("Blocks: %d, %d. 1\n", 1+(total-1)/BSIZE, NJBLOCK);
    dim3 nb(1+(total-1)/BSIZE, NJBLOCK,1);
    printf("Threads: %d, 1, 1\n", BSIZE);
    dim3 nt(BSIZE, 1 ,1);
    getchar();
    k_update_2d <<< nb, nt  >>> (d_move, d_new_a, d_new_j, d_p_r, d_p_v, d_m, n, total);
    cudaThreadSynchronize();
    std::cerr << "k_update_2d: " << cudaGetErrorString(cudaGetLastError()) <<  std::endl;

    cudaMemcpy(h_new_a,d_new_a,sizeof(double4) * n * NJBLOCK,cudaMemcpyDeviceToHost);
    cudaMemcpy(h_new_j,d_new_j,sizeof(double4) * n * NJBLOCK,cudaMemcpyDeviceToHost);



    for (int k = 0; k < n; k++) {
        int i = h_move[k];
        for (int j = 0; j < NJBLOCK; j++) {
            printf("%f %f %f\n",h_new_a[i + n * j].x, h_new_a[i + n * j].y, h_new_a[i + n*j].z );
        }
        getchar();
    }
}


/*
 * @fn gpu_energy()
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the energy of the system.
 */
__host__ double gpu_energy()
{
    k_energy <<< nblocks, nthreads >>> (d_p_r, d_p_v, d_ekin, d_epot, d_m, n);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_energy: " << cudaGetErrorString(cudaGetLastError()) << " | TpB: " << nthreads << " | BpG: " << nblocks << std::endl;
    #endif

    CUDA_SAFE_CALL(cudaMemcpy(h_ekin, d_ekin, f1_size,cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_epot, d_epot, f1_size,cudaMemcpyDeviceToHost));

    // Reduction on CPU
    ekin = 0.0f;
    epot = 0.0f;

    for (int i = 0; i < n; i++) {
        ekin += h_ekin[i];
        epot += h_epot[i];
    }
    return ekin + epot;
}

/*
 * @fn gpu_init_acc_jrk()
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the initial acceleration and jrk.
 */
__host__ void gpu_init_acc_jrk()
{
    int smem     = BSIZE * 2 * sizeof(double4);

    CUDA_SAFE_CALL(cudaMemcpy(d_r, h_r, d4_size, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v, h_v, d4_size, cudaMemcpyHostToDevice));

    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (d_r, d_v, d_a, d_j, d_m, n);
    //k_init_acc_jrk_simple <<< nblocks, nthreads >>> (d_r, d_v, d_a, d_j, d_m, n);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_init_acc_jrk: " << cudaGetErrorString(cudaGetLastError()) <<" | TpB: " << nthreads << " | BpG: " << nblocks << " | Smem: " << smem << std::endl;
    #endif

    CUDA_SAFE_CALL(cudaMemcpy(h_a, d_a, d4_size, cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_j, d_j, d4_size, cudaMemcpyDeviceToHost));
}


__host__ void gpu_predicted_pos_vel(float ITIME)
{
    k_predicted_pos_vel<<< nblocks, nthreads >>> (d_r,     d_v,     d_a,     d_j,
                                                  d_p_r,   d_p_v,   d_t, ITIME, n);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_predicted_pos_vel: " << cudaGetErrorString(cudaGetLastError()) <<" | TpB: " << nthreads << " | BpG: " << nblocks <<  std::endl;
    #endif

//    CUDA_SAFE_CALL(cudaMemcpy(h_p_r, d_p_r, d4_size, cudaMemcpyDeviceToHost));
//    CUDA_SAFE_CALL(cudaMemcpy(h_p_v, d_p_v, d4_size, cudaMemcpyDeviceToHost));

}

__host__ void gpu_update_acc_jrk_single(int total)
{
    int smem     = BSIZE * 2 * sizeof(double4);
    CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));

    for (int k = 0; k < total; k++)
    {
        int i = h_move[k];
        k_update_acc_jrk_single <<< nblocks, nthreads, smem  >>> (h_p_r[i], h_p_v[i], d_new_a, d_new_j, d_p_r, d_p_v, d_m, n, i);
        cudaThreadSynchronize();
        #ifdef KERNEL_ERROR_DEBUG
            std::cerr << "k_update_acc_jrk_single: " << cudaGetErrorString(cudaGetLastError()) <<" | TpB: " << nthreads << " | BpG: " << nblocks << " | Smem: " << smem <<  std::endl;
        #endif

        k_reduce<<< nblocks, nthreads, smem >>> (d_new_a, tmp_red, n);
        cudaThreadSynchronize();
        k_reduce<<< 1, nthreads, smem >>> (tmp_red, tmp_red + nblocks, nblocks);
        cudaMemcpy(h_a + i,tmp_red + nblocks,sizeof(double4),cudaMemcpyDeviceToHost);

        k_reduce<<< nblocks, nthreads, smem >>> (d_new_j, tmp_red, n);
        cudaThreadSynchronize();
        k_reduce<<< 1, nthreads, smem >>> (tmp_red, tmp_red + nblocks, nblocks);
        cudaMemcpy(h_j + i,tmp_red + nblocks,sizeof(double4),cudaMemcpyDeviceToHost);

        cudaMemcpy(h_a + i , d_a + i, sizeof(double4), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_j + i , d_j + i, sizeof(double4), cudaMemcpyDeviceToHost);
    }
}

__host__ void gpu_update_acc_jrk_simple(int total)
{
    CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));

    k_update_acc_jrk_simple <<< nblocks, nthreads >>> (d_p_r, d_p_v, d_a, d_j, d_m,
                                                            d_move, n, total);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_update_acc_jrk_simple: " << cudaGetErrorString(cudaGetLastError()) <<"TpB: " << nthreads << " BpG: " << nblocks << "Smem: " << smem <<  std::endl;
    #endif

    for (int i = 0; i < total; i++)
    {
        int j = h_move[i];
        CUDA_SAFE_CALL(cudaMemcpy(&h_a[j], &d_a[j], sizeof(double4), cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(&h_j[j], &d_j[j], sizeof(double4), cudaMemcpyDeviceToHost));
    }
}

/*
 * @fn gpu_update_acc_jrk()
 *
 * @param total amount of particles to update
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the new acceleration and jrk.
 */
__host__ void gpu_update_acc_jrk(int total)
{
    int smem     = BSIZE * 2 * sizeof(double4);

    CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));

    k_update_acc_jrk <<< nblocks, nthreads, smem >>> (d_p_r, d_p_v, d_a, d_j, d_m,
                                                            d_move, n, total);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_update_acc_jrk: " << cudaGetErrorString(cudaGetLastError()) <<"TpB: " << nthreads << " BpG: " << nblocks << "Smem: " << smem <<  std::endl;
    #endif

    for (int i = 0; i < total; i++)
    {
        int j = h_move[i];
        CUDA_SAFE_CALL(cudaMemcpy(&h_a[j], &d_a[j], sizeof(double4), cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(&h_j[j], &d_j[j], sizeof(double4), cudaMemcpyDeviceToHost));
    }
}
