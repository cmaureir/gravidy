#include "hermite.cuh"

void integrate_gpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    iterations   = 0;       // Iterations of the integration

    // Tmp setting nblocks and nthreads
    nthreads = BSIZE;
    nblocks = ceil(n/(float)nthreads);

    // Copying the input file information from the CPU to the GPU
    CUDA_SAFE_CALL(cudaMemcpy(d_r,  h_r,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v,  h_v,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_m,  h_m,  f1_size,cudaMemcpyHostToDevice));

    gpu_init_acc_jrk();   // Initial calculation of a and a1
    //init_acc_jrk();   // Initial calculation of a and a1

    // Copying a and a1 from the GPU to the CPU
    CUDA_SAFE_CALL(cudaMemcpy(h_f,  d_f,  sizeof(Forces) * n ,cudaMemcpyDeviceToHost));

    init_dt(&ATIME);  // Initial calculation of time-steps using simple equation
    //init_dt2(&ATIME); // Initial calculation of time-steps using complete equation

    energy_ini = gpu_energy(); // Initial calculation of the energy of the system
    energy_tmp = energy_ini;  // Saving initial energy, to calculate errors

    get_energy_log(ITIME, iterations, nsteps, out, energy_tmp); // First log of the integration

    float tmp_time = 0.0f;
    gpu_time = 0.0f;

    while (ITIME < int_time)
    {
        ITIME = ATIME;                         // New integration time
        nact = find_particles_to_move(ITIME);  // Find particles to move (nact)
        save_old(nact);                        // Save old information

        if (nact < n * alpha)
        {
            predicted_pos_vel(ITIME);
            update_acc_jrk(nact);
            correction_pos_vel(ITIME, nact);       // Correct r and v of nact particles
            cpu_iterations++;
        }
        else
        {
            tmp_time = (float)clock()/CLOCKS_PER_SEC;
            //CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));
            //CUDA_SAFE_CALL(cudaMemcpy(d_t,  h_t,  d1_size, cudaMemcpyHostToDevice));
            //gpu_predicted_pos_vel(ITIME);
            predicted_pos_vel(ITIME);
            CUDA_SAFE_CALL(cudaMemcpy(d_p, h_p, sizeof(Predictor) * n,cudaMemcpyHostToDevice));
            gpu_update(nact);     // Update a and a1 of nact particles
            //CUDA_SAFE_CALL(cudaMemcpy(h_f,  d_f,  sizeof(Forces) * n ,cudaMemcpyDeviceToHost));
            correction_pos_vel(ITIME, nact);       // Correct r and v of nact particles
            CUDA_SAFE_CALL(cudaMemcpy(d_r, h_r, d4_size, cudaMemcpyHostToDevice));
            CUDA_SAFE_CALL(cudaMemcpy(d_v, h_v, d4_size, cudaMemcpyHostToDevice));
            gpu_time += (float)clock()/CLOCKS_PER_SEC - tmp_time;
            gpu_iterations++;
        }

        next_itime(&ATIME);                    // Find next integration time

        if(std::ceil(ITIME) == ITIME)          // Print log in every integer ITIME
        {
           get_energy_log(ITIME, iterations, nsteps, out, gpu_energy());
        }

        nsteps += nact;                        // Update nsteps with nact
        iterations++;                          // Increase iterations
    }
}
