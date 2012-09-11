#include "dynamics_gpu.cuh"

struct sum_double4
{
    __host__ __device__
    double4 operator() (const double4 l, const double4 r)
    {
        return make_double4(l.x + r.x, l.y + r.y, l.z + r.z, 0.0);
    }
};

/*
 * @fn gpu_energy()
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the energy of the system.
 */
__host__
double
gpu_energy(bool type)
{
    int d4_size = sizeof(double4) * n;
    int f1_size = sizeof(float)  * n;

    int nthreads = BSIZE;
    int nblocks  = std::ceil(n/(float)nthreads);

    if(type)
    {
        cudaMemcpy(d_p_r, h_p_r, d4_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_p_v, h_p_v, d4_size, cudaMemcpyHostToDevice);
        k_energy <<< nblocks, nthreads >>> (d_p_r, d_p_v, d_ekin, d_epot, d_m, n);
    }
    else if(!type)
    {
        cudaMemcpy(d_r, h_r, d4_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_v, h_v, d4_size, cudaMemcpyHostToDevice);
        k_energy <<< nblocks, nthreads >>> (d_r, d_v, d_ekin, d_epot, d_m, n);
    }

    cudaThreadSynchronize();
//    std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;

    cudaMemcpy(h_ekin, d_ekin, f1_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_epot, d_epot, f1_size, cudaMemcpyDeviceToHost);

    ekin = 0;
    epot = 0;
    for (int i = 0; i < n; i++)
    {
        ekin += h_ekin[i];
        epot += h_epot[i];
    }

    return ekin + epot;

}

/*
 * @fn gpu_init_acc_jerk()
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the initial acceleration and jerk.
 */
__host__
void
gpu_init_acc_jerk()
{
    int d4_size = sizeof(double4) * n;
    int nthreads = BSIZE;
    int nblocks  = std::ceil(n/(float)nthreads);
    int smem     = BSIZE * 2 * sizeof(double4);

    cudaMemcpy(d_r, h_r, d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, h_v, d4_size, cudaMemcpyHostToDevice);

    //k_init_acc_jerk <<< nblocks, nthreads >>> (d_r, d_v, d_a, d_j, d_m, n);
    k_init_acc_jerk_tile <<< nblocks, nthreads, smem >>> (d_r, d_v, d_a, d_j, d_m, n);
    cudaThreadSynchronize();

    cudaMemcpy(h_a, d_a, d4_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_j, d_j, d4_size, cudaMemcpyDeviceToHost);
}

/*
 * @fn gpu_update_acc_jerk()
 *
 * @param total amount of particles to update
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the new acceleration and jerk.
 */
__host__
void
gpu_update_acc_jerk_simple(int total)
{
    int d4_size = sizeof(double4) * n;
    int i1_size = sizeof(int) * n;
    int nthreads = BSIZE;
    int nblocks  = std::ceil(n/(float)nthreads);

    cudaMemcpy(d_p_r,  h_p_r,  d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p_v,  h_p_v,  d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice);

    k_update_acc_jerk_simple <<< nblocks, nthreads      >>> (d_p_r, d_p_v, d_a, d_j, d_m,
                                                             d_move, n, total);
    cudaThreadSynchronize();

    cudaMemcpy(h_a, d_a, d4_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_j, d_j, d4_size, cudaMemcpyDeviceToHost);
}

__host__
void
gpu_update_acc_jerk_single(int total)
{
    int d4_size = sizeof(double4) * n;
    int nthreads = BSIZE;
    int nblocks  = std::ceil(n/(float)nthreads);
    int smem     = BSIZE * 2 * sizeof(double4);
    double4 new_a = {0.0,0.0,0.0,0.0};
    double4 new_j = {0.0,0.0,0.0,0.0};
    thrust::device_ptr<double4> dptr_a(d_new_a);
    thrust::device_ptr<double4> dptr_j(d_new_j);

    cudaMemcpy(d_p_r,    h_p_r,    d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p_v,    h_p_v,    d4_size, cudaMemcpyHostToDevice);

    for (int k = 0; k < total; k++) {
        int i = h_move[k];
        k_update_acc_jerk_single <<< nblocks, nthreads, smem  >>> (h_p_r[i], h_p_v[i], d_new_a, d_new_j, d_p_r, d_p_v, d_m, n, i);
        
        cudaThreadSynchronize();

        try
        {
          new_a = thrust::reduce(dptr_a, dptr_a + n, make_double4(0.0,0.0,0.0,0.0), sum_double4());
          new_j = thrust::reduce(dptr_j, dptr_j + n, make_double4(0.0,0.0,0.0,0.0), sum_double4());
        }
        catch(std::bad_alloc &e)
        {
            std::cerr << "Ran out of memory while sorting" << std::endl;
            exit(-1);
        }
        catch(thrust::system_error &e)
        {
            std::cerr << "Some other error happened during sort: " << e.what() << std::endl;
            exit(-1);
        }

        h_a[i] = new_a;
        h_j[i] = new_j;
    }
}

__host__
void
gpu_update_acc_jerk_tile(int total)
{
    int d4_size = sizeof(double4) * n;
    int i1_size = sizeof(int) * n;
    int nthreads = BSIZE;
    int nblocks  = std::ceil(n/(float)nthreads);
    int smem     = BSIZE * 2 * sizeof(double4);

    cudaMemcpy(d_p_r,  h_p_r,  d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p_v,  h_p_v,  d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice);

    k_update_acc_jerk_tile <<< nblocks, nthreads, smem >>> (d_p_r, d_p_v, d_a, d_j, d_m,
                                                            d_move, n, total);
    cudaThreadSynchronize();

    cudaMemcpy(h_a, d_a, d4_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_j, d_j, d4_size, cudaMemcpyDeviceToHost);
}

/*
 * @fn gpu_update_acc_jerk()
 *
 * @param total amount of particles to update
 *
 * @brief
 *  Host function which call the kernel to calculate
 *  the new acceleration and jerk.
 */
__host__
void
gpu_correction_pos_vel(double ITIME, int total)
{
    int d4_size = sizeof(double4) * n;
    int f1_size = sizeof(double) * n;
    int i1_size = sizeof(int) * total;
    // TO DO
    // Fix to work only with the particles who need to be moved
    int nthreads = BSIZE;
    int nblocks  = std::ceil(n/(float)nthreads);

    cudaMemcpy(d_r,     h_r,     d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_v,     h_v,     d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_a,     h_a,     d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_j,     h_j,     d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_old_a, h_old_a, d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_old_j, h_old_j, d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p_r,   h_p_r,   d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p_v,   h_p_v,   d4_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_t,     h_t,     f1_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dt,    h_dt,    f1_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_move,  h_move,  i1_size, cudaMemcpyHostToDevice);

    k_correction_pos_vel<<< nblocks, nthreads >>> (d_r,     d_v,     d_a,     d_j,
                                                   d_old_a, d_old_j, 
                                                   d_p_r,   d_p_v,   d_t,     d_dt,
                                                   ITIME,   d_move,  total);
    cudaThreadSynchronize();

    cudaMemcpy(h_r,  d_r,  d4_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_v,  d_v,  d4_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_t,  d_t,  f1_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_dt, d_dt, f1_size, cudaMemcpyDeviceToHost);
}
