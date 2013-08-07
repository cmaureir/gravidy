#include "hermite.cuh"

void integrate_gpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    iterations   = 0;       // Iterations of the integration
    static long long interactions = 0;

    // General setting nblocks and nthreads
    nthreads = BSIZE;
    nblocks = ceil(n/(float)nthreads);

    // Copying the input file information from the CPU to the GPU
    CUDA_SAFE_CALL(cudaMemcpy(d_r,  h_r,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v,  h_v,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_m,  h_m,  f1_size,cudaMemcpyHostToDevice));

    // Initial calculation of a and a1 on the GPU
    gpu_init_acc_jrk();

    // Copying a and a1 from the GPU to the CPU
    CUDA_SAFE_CALL(cudaMemcpy(h_f,  d_f,  sizeof(Forces) * n ,cudaMemcpyDeviceToHost));

    // Initial calculation of time-steps using simple equation
    init_dt(&ATIME);

    // Initial calculation of the energy of the system
    energy_ini = gpu_energy();

    // Saving initial energy, to calculate errors
    energy_tmp = energy_ini;

    // Print First log
    get_energy_log(ITIME, iterations, interactions, nsteps, out, energy_tmp);


    while (ITIME < itime)
    {
        // New integration time
        ITIME = ATIME;

        // Find particles to update this iteration (nact)
        nact = find_particles_to_move(ITIME);

        // Save current values of the a and a1 to perform correction later
        save_old(nact);

        // Prediction of the nact particles
        predicted_pos_vel(ITIME);

        // Update of the gravitational interactions of the nact particles
        gpu_update(nact);

        // Correction of the nact particles
        correction_pos_vel(ITIME, nact);


        // Update the amount of interactions counter
        interactions += nact * n;

        // Find the next integration time
        next_itime(&ATIME);

        // Print log every integer ITIME
        if(std::ceil(ITIME) == ITIME)
        //if(nact == n)          // Print log in every integer ITIME
        {
           get_energy_log(ITIME, iterations, interactions, nsteps, out, gpu_energy());
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        iterations++;
    }
}
