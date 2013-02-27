#include "hermite.hpp"

void integrate_cpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    iterations   = 0;       // Iterations of the integration

    init_acc_jrk();   // Initial calculation of a and a1
    init_dt(&ATIME);  // Initial calculation of time-steps using simple equation
    //init_dt2(&ATIME); // Initial calculation of time-steps using complete equation

    energy_ini = energy();   // Initial calculation of the energy of the system
    energy_tmp = energy_ini; // Saving initial energy, to calculate errors

    get_energy_log(ITIME, iterations, nsteps, out, energy_tmp); // First log of the integration

    while (ITIME < int_time)
    {
        ITIME = ATIME;                         // New integration time
        nact = find_particles_to_move(ITIME);  // Find particles to move (nact)
        save_old(nact);                        // Save old information

        predicted_pos_vel(ITIME);              // Predict all the particles
        update_acc_jrk(nact);                  // Update a and a1 of nact particles
        correction_pos_vel(ITIME, nact);       // Correct r and v of nact particles

        next_itime(&ATIME);                    // Find next integration time

        if(std::ceil(ITIME) == ITIME)          // Print log in every integer ITIME
        {
           get_energy_log(ITIME, iterations, nsteps, out, energy());
           //print_all(n,ITIME);
        }

        nsteps += nact;                        // Update nsteps with nact
        iterations++;                          // Increase iterations
        print_all(n,ITIME);
        if(iterations == 10) ITIME = int_time;
    }
}
