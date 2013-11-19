#include "../NbodySystem.hpp"

double NbodySystem::get_energy_gpu()
{
    CUDA_SAFE_CALL(cudaMemcpy(d_r,  h_r,  sizeof(double4) * n,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v,  h_v,  sizeof(double4) * n,cudaMemcpyHostToDevice));

    //gpu_timer_start();
    int nthreads = BSIZE;
    int nblocks = std::ceil(n/(float)nthreads);
    k_energy <<< nblocks, nthreads >>> (d_r, d_v, d_ekin, d_epot, n, e2);
    cudaThreadSynchronize();
    //float msec = gpu_timer_stop("k_energy");
    //get_kernel_error();

    CUDA_SAFE_CALL(cudaMemcpy(h_ekin, d_ekin, sizeof(double) * n,cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_epot, d_epot, sizeof(double) * n,cudaMemcpyDeviceToHost));

    // Reduction on CPU
    en.kinetic = 0.0;
    en.potential = 0.0;

    for (int i = 0; i < n; i++) {
        en.kinetic   += h_ekin[i];
        en.potential += h_epot[i];
    }
    return en.kinetic + en.potential;
}

__global__ void k_energy(double4 *r,
                         double4 *v,
                         double *ekin,
                         double *epot,
                         int n,
                         double e2)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j;
    double ekin_tmp = 0.0;

    if (i < n)
    {
        double epot_tmp = 0.0;
        for (j = i+1; j < n; j++)
        {
            double rx = r[j].x - r[i].x;
            double ry = r[j].y - r[i].y;
            double rz = r[j].z - r[i].z;
            double r2 = rx*rx + ry*ry + rz*rz + e2;

            epot_tmp -= (r[i].w * r[j].w) * rsqrt(r2);
        }

        double vx = v[i].x * v[i].x;
        double vy = v[i].y * v[i].y;
        double vz = v[i].z * v[i].z;
        double v2 = vx + vy + vz;

        ekin_tmp = 0.5 * r[i].w * v2;

        ekin[i] = ekin_tmp;
        epot[i] = epot_tmp;
    }
}
