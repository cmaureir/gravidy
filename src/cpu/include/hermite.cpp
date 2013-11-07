#include "hermite.hpp"

void integrate_cpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    iterations   = 0;       // Iterations of the integration
    static long long interactions = 0;

    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    init_acc_jrk();   // Initial calculation of a and a1
    init_dt(&ATIME);  // Initial calculation of time-steps using simple equation

    energy_ini = energy();   // Initial calculation of the energy of the system
    energy_tmp = energy_ini; // Saving initial energy, to calculate errors

    get_energy_log(ITIME, iterations, interactions, nsteps, out, energy_tmp);

    while (ITIME < itime)
    {
        ITIME = ATIME;                         // New integration time
        nact = find_particles_to_move(ITIME);  // Find particles to move (nact)
        save_old(nact);                        // Save old information

        predicted_pos_vel(ITIME);              // Predict all the particles
        update_acc_jrk(nact);                  // Update a and a1 of nact particles
        correction_pos_vel(ITIME, nact);       // Correct r and v of nact particles


        // Update the amount of interactions counter
        interactions += nact * n;

        // Find the next integration time
        next_itime(&ATIME);

        // Print log every integer ITIME
        if(std::ceil(ITIME) == ITIME)
        //if(nact == n)          // Print log in every integer ITIME
        {
           get_energy_log(ITIME, iterations, interactions, nsteps, out, energy());
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        iterations++;
    }
}
