#include "dynamics_gpu.cuh"

__host__ void gpu_init_acc_jrk()
{
    int smem = BSIZE * 2 * sizeof(double4);
    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (d_r, d_v, d_a, d_a1, d_m, n);
    cudaThreadSynchronize();
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
            omp_get_wtime() - ini_time,
            energy_end,
            relative_error,
            energy_tmp/(int)ITIME);
}

__host__ double gpu_energy()
{
    k_energy <<< nblocks, nthreads >>> (d_r, d_v, d_ekin, d_epot, d_m, n);
    cudaThreadSynchronize();
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

__host__ void gpu_predicted_pos_vel(float ITIME)
{
    k_predicted_pos_vel<<< nblocks, nthreads >>> (d_r, d_v, d_a, d_a1, d_p_r, d_p_v,
                                                  d_t, ITIME, n);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_predicted_pos_vel: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}


__host__ void gpu_update_acc_jrk_simple(int total)
{
    //CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));

    int smem = BSIZE * 2 * sizeof(double4);
    k_update_acc_jrk_simple <<< nblocks, nthreads, smem >>> (d_p_r, d_p_v, d_a, d_a1,
                                                             d_m, d_move, n, total);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_update_acc_jrk_simple: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}


//__host__ void gpu_update_acc_jrk(int total)
//{
//    //int smem     = BSIZE * 2 * sizeof(double4);
//    CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));
//
//    for (int k = 0; k < total; k++)
//    {
//        int i = h_move[k];
//        k_update_acc_jrk_single <<< nblocks, nthreads  >>> (d_new_a, d_new_j,
//                                                           d_p_r, d_p_v, d_m, n, i);
//
//        cudaThreadSynchronize();
//        #ifdef KERNEL_ERROR_DEBUG
//            std::cerr << "k_update_acc_jrk_single: " << std::endl;
//            std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
//        #endif
//
//        cudaMemcpy(h_new_a,d_new_a,d4_size,cudaMemcpyDeviceToHost);
//        cudaMemcpy(h_new_j,d_new_j,d4_size,cudaMemcpyDeviceToHost);
//
//        h_a[i].x = 0.0;
//        h_a[i].y = 0.0;
//        h_a[i].z = 0.0;
//
//        h_a1[i].x = 0.0;
//        h_a1[i].y = 0.0;
//        h_a1[i].z = 0.0;
//
//
//        for (int iii = 0; iii < n; iii++) {
//            h_a[i] += h_new_a[iii];
//            h_a1[i] += h_new_j[iii];
//        }
//
////        k_reduce<<< nblocks, nthreads, smem >>> (d_new_a, tmp_red, n);
////        cudaThreadSynchronize();
////        k_reduce<<< 1, nthreads, smem >>> (tmp_red, tmp_red + nblocks, nblocks);
////        cudaMemcpy(h_a + i,tmp_red + nblocks,sizeof(double4),cudaMemcpyDeviceToHost);
////
////        k_reduce<<< nblocks, nthreads, smem >>> (d_new_j, tmp_red, n);
////        cudaThreadSynchronize();
////        k_reduce<<< 1, nthreads, smem >>> (tmp_red, tmp_red + nblocks, nblocks);
////        cudaMemcpy(h_a1 + i,tmp_red + nblocks,sizeof(double4),cudaMemcpyDeviceToHost);
////
//        //cudaMemcpy(h_a + i , d_a + i, sizeof(double4), cudaMemcpyDeviceToHost);
//        //cudaMemcpy(h_a1 + i , d_a1 + i, sizeof(double4), cudaMemcpyDeviceToHost);
//    }
//}
//__host__ void gpu_update_2d(int total)
//{
//    // i-threads = number of particles to move (total)
//    // j-threads = particles per block to calculate force on i (~ 32)
//    ////////ENVIAR SOLO PARTICULAS QUE SE MUEVEN
//    cudaMemcpy(d_move, h_move, n * sizeof(int), cudaMemcpyHostToDevice);
//
//    dim3 nb(1+(total-1)/BSIZE, NJBLOCK,1);
//    dim3 nt(BSIZE, 1 ,1);
//    printf("Moving: %4d using Blocks(%3d, %3d. 1)  Threads(%4d, 1, 1)\n",
//            total, nb.x, nb.y, nt.x);
//    k_update_2d <<< nb, nt  >>> (d_move, d_new_a, d_new_j, d_p_r, d_p_v, d_m, n, total);
//    cudaThreadSynchronize();
//    std::cerr << "k_update_2d: " << cudaGetErrorString(cudaGetLastError()) <<  std::endl;
//
//    cudaMemcpy(h_new_a,d_new_a,sizeof(double4) * n * NJBLOCK,cudaMemcpyDeviceToHost);
//    cudaMemcpy(h_new_j,d_new_j,sizeof(double4) * n * NJBLOCK,cudaMemcpyDeviceToHost);
//
//    for (int k = 0; k < n; k++) {
//        int i = h_move[k];
//        for (int j = 0; j < NJBLOCK; j++) {
//            printf("%f %f %f\n",h_new_a[i + n * j].x, h_new_a[i + n * j].y, h_new_a[i + n*j].z );
//        }
//        printf("\n");
//    }
//
//        getchar();
//}
//


/*
 * @fn gpu_update_acc_jrk()
 *
 * @param total amount of particles to update
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the new acceleration and jrk.
 */
 /*
__host__ void gpu_update_acc_jrk(int total)
{
    int smem     = BSIZE * 2 * sizeof(double4);

    CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));

    k_update_acc_jrk <<< nblocks, nthreads, smem >>> (d_p_r, d_p_v, d_a, d_a1, d_m,
                                                            d_move, n, total);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_update_acc_jrk: " << cudaGetErrorString(cudaGetLastError()) <<"TpB: " << nthreads << " BpG: " << nblocks << "Smem: " << smem <<  std::endl;
    #endif

    for (int i = 0; i < total; i++)
    {
        int j = h_move[i];
        CUDA_SAFE_CALL(cudaMemcpy(&h_a[j], &d_a[j], sizeof(double4), cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(&h_a1[j], &d_a1[j], sizeof(double4), cudaMemcpyDeviceToHost));
    }
}
*/
