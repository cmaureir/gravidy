#include "dynamics_gpu.cuh"

//struct sum_double4
//{
//    __host__ __device__
//    double4 operator() (const double4 l, const double4 r)
//    {
//        return make_double4(l.x + r.x, l.y + r.y, l.z + r.z, 0);
//    }
//};


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


__host__ void gpu_init_dt(float *ATIME)
{
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    for (int i = 0; i < n; i++) {
        float tmp_dt = ETA_S *
                 magnitude(h_a[i].x, h_a[i].y, h_a[i].z) /
                 magnitude(h_j[i].x, h_j[i].y, h_j[i].z);

        /* Adjusting to block timesteps */
        tmp_dt = pow(2,(int)((log(tmp_dt)/log(2.0))-1));
        if (tmp_dt < D_TIME_MIN)
            tmp_dt = D_TIME_MIN;
        else if (tmp_dt > D_TIME_MAX)
            tmp_dt = D_TIME_MAX;

        h_dt[i] = tmp_dt;
        h_t[i] = 0.0;

        // Obtaining the first integration time
        if(tmp_dt < *ATIME)
        {
            *ATIME = tmp_dt;
        }
    }
    CUDA_SAFE_CALL(cudaMemcpy(d_dt,h_dt,f1_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_t, h_t, f1_size,cudaMemcpyHostToDevice));
    // TMP, first send
    CUDA_SAFE_CALL(cudaMemcpy(d_p_r, h_p_r, d4_size, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_p_v, h_p_v, d4_size, cudaMemcpyHostToDevice));

}


__host__ void gpu_next_itime(float *ATIME)
{
    *ATIME = 1.0e10;
    for (int i = 0; i < n; i++)
    {
        float time = h_t[i] + h_dt[i];
        if(time < *ATIME)
            *ATIME = time;
    }
}

__host__ int gpu_find_particles_to_move(float ITIME)
{
    int j = 0;
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        h_move[i] = -1;
        if (h_t[i] + h_dt[i] == ITIME)
        {
            h_move[j] = i;
            j++;
        }
    }
    return j;
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

    CUDA_SAFE_CALL(cudaMemcpy(h_p_r, d_p_r, d4_size, cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_p_v, d_p_v, d4_size, cudaMemcpyDeviceToHost));

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

__host__
void
gpu_correction_pos_vel(float ITIME, int total)
{
//    // TO DO
//    // Fix to work only with the particles who need to be moved
//    cudaMemcpy(d_dt,    h_dt,    f1_size, cudaMemcpyHostToDevice);
//    k_correction_pos_vel<<< std::ceil(total/(float)32), 32 >>> (d_r, d_v, d_a, d_j,
//                                                                d_old_a,  d_old_j,
//                                                                d_p_r,    d_p_v,
//                                                                d_t,     d_dt,
//                                                   ITIME,   d_move,  total);
//    cudaThreadSynchronize();
//    #ifdef KERNEL_ERROR_DEBUG
//        std::cerr << "k_correction_pos_vel: " << cudaGetErrorString(cudaGetLastError()) <<"TpB: " << nthreads << " BpG: " << nblocks <<  std::endl;
//    #endif
//
//    //cudaMemcpy(h_r,  d_r,  d4_size, cudaMemcpyDeviceToHost);
//    //cudaMemcpy(h_v,  d_v,  d4_size, cudaMemcpyDeviceToHost);
//    cudaMemcpy(h_t,  d_t,  f1_size, cudaMemcpyDeviceToHost);
//    cudaMemcpy(h_dt, d_dt, f1_size, cudaMemcpyDeviceToHost);
//
// on CPU


    //cudaMemcpy(h_r  , d_r ,  d4_size,cudaMemcpyDeviceToHost);
    //cudaMemcpy(h_v  , d_v ,  d4_size,cudaMemcpyDeviceToHost);
    //cudaMemcpy(h_p_r, d_p_r ,d4_size,cudaMemcpyDeviceToHost);
    //cudaMemcpy(h_p_v, d_p_v ,d4_size,cudaMemcpyDeviceToHost);


    for (int k = 0; k < total; k++)
    {
        int i = h_move[k];

        cudaMemcpy(h_r   + i, d_r + i ,  sizeof(double4),cudaMemcpyDeviceToHost);
        cudaMemcpy(h_v   + i, d_v + i ,  sizeof(double4),cudaMemcpyDeviceToHost);
        cudaMemcpy(h_p_r + i, d_p_r + i, sizeof(double4),cudaMemcpyDeviceToHost);
        cudaMemcpy(h_p_v + i, d_p_v + i, sizeof(double4),cudaMemcpyDeviceToHost);

        float dt1 = ITIME - h_t[i];
        h_t[i] = ITIME;
        float dt2 = dt1 * dt1;
        float dt3 = dt2 * dt1;
        float dt4 = dt2 * dt2;
        float dt5 = dt4 * dt1;

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
        double abs_a1   = sqrt((h_a[i].x * h_a[i].x) + (h_a[i].y * h_a[i].y) + (h_a[i].z * h_a[i].z));
        // |j_{1,i}|
        double abs_j1   = sqrt((h_j[i].x * h_j[i].x) + (h_j[i].y * h_j[i].y) + (h_j[i].z * h_j[i].z));
        // |j_{1,i}|^{2}
        double abs_j12  = abs_j1 * abs_j1;
        // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
        double abs_a1_3 = sqrt((ax0_3 * ax0_3) + (ay0_3 * ay0_3) + (az0_3 * az0_3));
        // |a_{1,i}^{(2)}|
        double abs_a1_2 = sqrt((ax1_2 * ax1_2) + (ay1_2 * ay1_2) + (az1_2 * az1_2));
        // |a_{1,i}^{(2)}|^{2}
        double abs_a1_22  = abs_a1_2 * abs_a1_2;

        float tmp_dt = sqrt(ETA_N * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));

        /* Adjusting to block timesteps */
        if (tmp_dt < h_dt[i])
            tmp_dt = h_dt[i]/2;
        else if (tmp_dt > 2 * h_dt[i])
            tmp_dt = 2 * h_dt[i];
        else
            tmp_dt = h_dt[i];

        if (tmp_dt < D_TIME_MIN)
            tmp_dt = D_TIME_MIN;
        else if (tmp_dt > D_TIME_MAX)
            tmp_dt = D_TIME_MAX;

        h_dt[i] = tmp_dt;

        cudaMemcpy(d_r  + i , h_r  + i , sizeof(double4), cudaMemcpyHostToDevice);
        cudaMemcpy(d_v  + i , h_v  + i , sizeof(double4), cudaMemcpyHostToDevice);
        cudaMemcpy(d_t  + i , h_t  + i , sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dt + i,  h_dt + i , sizeof(float), cudaMemcpyHostToDevice);
    }
//        cudaMemcpy(d_r  , h_r , d4_size, cudaMemcpyHostToDevice);
//        cudaMemcpy(d_v  , h_v , d4_size, cudaMemcpyHostToDevice);
//        cudaMemcpy(d_t  , h_t , f1_size, cudaMemcpyHostToDevice);
//        cudaMemcpy(d_dt,  h_dt, f1_size, cudaMemcpyHostToDevice);



}
