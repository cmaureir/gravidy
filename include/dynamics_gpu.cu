#include "dynamics_gpu.cuh"

void gpu_correction_pos_vel(double ITIME, int total)
{

    for (int k = 0; k < total; k++)
    {
        int i = h_move[k];

        CUDA_SAFE_CALL(cudaMemcpy(h_p_r + i,  d_p_r + i , sizeof(double4),cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(h_p_v + i,  d_p_v + i , sizeof(double4),cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(h_a   + i,  d_a   + i , sizeof(double4),cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(h_j   + i,  d_j   + i , sizeof(double4),cudaMemcpyDeviceToHost));


        double dt1 = h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        double ax0_2 = (-6 * (h_old_a[i].x - h_a[i].x ) - dt1 * (4 * h_old_j[i].x + 2 * h_j[i].x) ) / dt2;
        double ay0_2 = (-6 * (h_old_a[i].y - h_a[i].y ) - dt1 * (4 * h_old_j[i].y + 2 * h_j[i].y) ) / dt2;
        double az0_2 = (-6 * (h_old_a[i].z - h_a[i].z ) - dt1 * (4 * h_old_j[i].z + 2 * h_j[i].z) ) / dt2;

        // Acceleration 3rd derivate
        double ax0_3 = (12 * (h_old_a[i].x - h_a[i].x ) + 6 * dt1 * (h_old_j[i].x + h_j[i].x) ) / dt3;
        double ay0_3 = (12 * (h_old_a[i].y - h_a[i].y ) + 6 * dt1 * (h_old_j[i].y + h_j[i].y) ) / dt3;
        double az0_3 = (12 * (h_old_a[i].z - h_a[i].z ) + 6 * dt1 * (h_old_j[i].z + h_j[i].z) ) / dt3;

        // Correcting position
        h_r[i].x = h_p_r[i].x + (dt4/24)*ax0_2 + (dt5/120)*ax0_3;
        h_r[i].y = h_p_r[i].y + (dt4/24)*ay0_2 + (dt5/120)*ay0_3;
        h_r[i].z = h_p_r[i].z + (dt4/24)*az0_2 + (dt5/120)*az0_3;

        // Correcting velocity
        h_v[i].x = h_p_v[i].x + (dt3/6)*ax0_2 + (dt4/24)*ax0_3;
        h_v[i].y = h_p_v[i].y + (dt3/6)*ay0_2 + (dt4/24)*ay0_3;
        h_v[i].z = h_p_v[i].z + (dt3/6)*az0_2 + (dt4/24)*az0_3;

        // Timestep update

        // Calculating a_{1,i}^{(2)} = a_{0,i}^{(2)} + dt * a_{0,i}^{(3)}
        double ax1_2 = ax0_2 + dt1 * ax0_3;
        double ay1_2 = ay0_2 + dt1 * ay0_3;
        double az1_2 = az0_2 + dt1 * az0_3;

        // |a_{1,i}|
        double abs_a1 = get_magnitude(h_a[i].x, h_a[i].y, h_a[i].z);
        // |j_{1,i}|
        double abs_j1 = get_magnitude(h_j[i].x, h_j[i].y, h_j[i].z);
        // |j_{1,i}|^{2}
        double abs_j12  = abs_j1 * abs_j1;
        // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
        double abs_a1_3 = get_magnitude(ax0_3, ay0_3, az0_3);
        // |a_{1,i}^{(2)}|
        double abs_a1_2 = get_magnitude(ax1_2, ay1_2, az1_2);
        // |a_{1,i}^{(2)}|^{2}
        double abs_a1_22  = abs_a1_2 * abs_a1_2;

        double tmp_dt = sqrt(ETA_N * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));

        h_t[i] = ITIME;
        tmp_dt = normalize_dt(tmp_dt, h_dt[i], h_t[i], i);
        h_dt[i] = tmp_dt;

        CUDA_SAFE_CALL(cudaMemcpy(d_t + i,  h_t + i, sizeof(double),cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_r + i,  h_r + i, sizeof(double4),cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_v + i,  h_v + i, sizeof(double4),cudaMemcpyHostToDevice));
    }
}
__host__ void gpu_send_initial_data()
{
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

__host__ void gpu_recv_initial_data()
{
}

__host__ void gpu_init_acc_jrk()
{
    int smem = BSIZE * 2 * sizeof(double4);
    k_init_acc_jrk <<< nblocks, nthreads, smem >>> (d_r, d_v, d_a, d_j, d_m, n);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_init_acc_jrk: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
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
    printf("%.10f %.10f\n", ekin, epot);
    return ekin + epot;
}

__host__ void gpu_predicted_pos_vel(float ITIME)
{


    k_predicted_pos_vel<<< nblocks, nthreads >>> (d_r, d_v, d_a, d_j, d_p_r, d_p_v,
                                                  d_t, ITIME, n);
    cudaThreadSynchronize();
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "k_predicted_pos_vel: " << std::endl;
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}


__host__ void gpu_update_acc_jrk_simple(int total)
{
    CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));

    int smem = BSIZE * 2 * sizeof(double4);
    k_update_acc_jrk_simple <<< nblocks, nthreads, smem >>> (d_p_r, d_p_v, d_a, d_j,
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
//        h_j[i].x = 0.0;
//        h_j[i].y = 0.0;
//        h_j[i].z = 0.0;
//
//
//        for (int iii = 0; iii < n; iii++) {
//            h_a[i] += h_new_a[iii];
//            h_j[i] += h_new_j[iii];
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
////        cudaMemcpy(h_j + i,tmp_red + nblocks,sizeof(double4),cudaMemcpyDeviceToHost);
////
//        //cudaMemcpy(h_a + i , d_a + i, sizeof(double4), cudaMemcpyDeviceToHost);
//        //cudaMemcpy(h_j + i , d_j + i, sizeof(double4), cudaMemcpyDeviceToHost);
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
*/

double gpu_normalize_dt(double new_dt, double old_dt, double t, int i)
{
    if (new_dt <= old_dt/8)
    {
        new_dt = D_TIME_MIN;
    }
    else if ( old_dt/8 < new_dt && new_dt <= old_dt/4)
    {
        new_dt = old_dt / 8;
    }
    else if ( old_dt/4 < new_dt && new_dt <= old_dt/2)
    {
        new_dt = old_dt / 4;
    }
    else if ( old_dt/2 < new_dt && new_dt <= old_dt)
    {
        new_dt = old_dt / 2;
    }
    else if ( old_dt < new_dt && new_dt <= old_dt * 2)
    {
        new_dt = old_dt;
    }
    else if (2 * old_dt < new_dt)
    {
        float val = t/(2 *old_dt);
        if(std::ceil(val) == val)
        {
            new_dt = 2 * old_dt;
        }
        else
        {
            new_dt = old_dt;
        }
    }
    else
    {
        fprintf(stderr, "gravidy: Undefined treatment for the time-step of (%d)",i);
        getchar();
    }

    if (new_dt < D_TIME_MIN)
    {
        new_dt = D_TIME_MIN;
    }
    else if (new_dt > D_TIME_MAX)
    {
        new_dt = D_TIME_MAX;
    }

    return new_dt;
}
