#include "hermite.cuh"

void integrate_gpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    iterations   = 0;       // Iterations of the integration

    // Copying the input file information from the CPU to the GPU
    CUDA_SAFE_CALL(cudaMemcpy(d_r,  h_r,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v,  h_v,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_m,  h_m,  f1_size,cudaMemcpyHostToDevice));

    gpu_init_acc_jrk();   // Initial calculation of a and a1
    //init_acc_jrk();   // Initial calculation of a and a1

    // Copying a and a1 from the GPU to the CPU
    CUDA_SAFE_CALL(cudaMemcpy(h_a,  d_a,  d4_size,cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_a1, d_a1, d4_size,cudaMemcpyDeviceToHost));

    init_dt(&ATIME);  // Initial calculation of time-steps using simple equation
    //init_dt2(&ATIME); // Initial calculation of time-steps using complete equation

    energy_ini = gpu_energy(); // Initial calculation of the energy of the system
    energy_tmp = energy_ini;  // Saving initial energy, to calculate errors

    get_energy_log(ITIME, iterations, nsteps, out); // First log of the integration

    CUDA_SAFE_CALL(cudaMemcpy(d_t,  h_t,  d1_size,cudaMemcpyHostToDevice));
    while (ITIME < int_time)
    {
        ITIME = ATIME;                         // New integration time
        nact = find_particles_to_move(ITIME);  // Find particles to move (nact)
        save_old(nact);                        // Save old information

        // Copying the index of the particles to move from the CPU to the GPU
        CUDA_SAFE_CALL(cudaMemcpy(d_move, h_move, i1_size, cudaMemcpyHostToDevice));


        gpu_predicted_pos_vel(ITIME);          // Predict all the particles
        gpu_update_acc_jrk_simple(nact);     // Update a and a1 of nact particles
        //update_acc_jrk(nact);
        //getchar();

        CUDA_SAFE_CALL(cudaMemcpy(h_p_r, d_p_r, d4_size,cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(h_p_v, d_p_v, d4_size,cudaMemcpyDeviceToHost));
        //// Copying a and a1 from the GPU to the CPU
        CUDA_SAFE_CALL(cudaMemcpy(h_a,  d_a,  d4_size,cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(h_a1, d_a1, d4_size,cudaMemcpyDeviceToHost));

        correction_pos_vel(ITIME, nact);       // Correct r and v of nact particles

        // Copying r, v and t from CPU to GPU
        CUDA_SAFE_CALL(cudaMemcpy(d_r, h_r, d4_size, cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_v, h_v, d4_size, cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_t, h_t, d1_size, cudaMemcpyHostToDevice));

        next_itime(&ATIME);                    // Find next integration time

        if(std::ceil(ITIME) == ITIME)          // Print log in every integer ITIME
        {
           get_energy_log(ITIME, iterations, nsteps, out);
        }

        //printf("%f\n", ITIME);
        nsteps += nact;                        // Update nsteps with nact
        iterations++;                          // Increase iterations
//        printf("# %d %d %.10f\n", iterations, nact, ITIME);
//        print_all(n, ITIME);
//       if(iterations == 50) ITIME = int_time;
    }
    print_all(n,0);
}
