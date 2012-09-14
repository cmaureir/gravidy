#include "hermite.hpp"
#define ITE_MAX (1e4)
#define OUT (1e3)

void
integrate_gpu()
{

    float ATIME = 1.0e+10; // Actual integration time
    float ITIME = 0.0;     // Integration time
    int total;

    energy_total = 0.0f;

    // Initial Acceleration and Jerk calculation on CPU
    init_time = omp_get_wtime();
    gpu_init_acc_jerk();
    init_time = omp_get_wtime() - init_time;

    // Initialazing block time steps
    init_dt(&ATIME);

    energy_ini = gpu_energy(0);      // Get initial energy
    //energy_ini = initial_energy();      // Get initial energy
    fprintf(stderr, "Initial energy: %.10f\n", energy_ini);
    iterations = 0;

    while (ITIME < int_time && iterations < ITE_MAX)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old();
        predicted_pos_vel(ITIME);
        if(total < BSIZE * 4)
        {
            gpu_update_acc_jerk_single(total);
        }
        else
        {
            //gpu_update_acc_jerk_simple(total);
            gpu_update_acc_jerk_tile(total);
        }
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;

        if(iterations % (int)OUT == 0)
        {
            //energy_end = gpu_energy(1);
            //fprintf(stderr, "%.10f %.10f %.10f %.5e\n", ITIME, energy_end, energy_end - energy_ini, (energy_end - energy_ini)/energy_ini);
            //print_positions(100);
            std::cout << ITIME << " " << iterations << std::endl;
        }

    }
    //energy_end = gpu_energy(1);
    energy_end = energy();
    fprintf(stderr, "Final energy: %.10f\n", energy_end);
    fprintf(stderr, "DeltaE: %.10e\n", (energy_end - energy_ini)/energy_ini);
}

void
integrate_cpu()
{
    float ATIME = 1.0e+10; // Actual integration time
    float ITIME = 0.0;     // Integration time
    int total = 0;

    energy_total = 0.0f;

    init_time = omp_get_wtime();
    init_acc_jerk();
    init_time = omp_get_wtime() - init_time;

    init_dt(&ATIME);
    energy_ini = initial_energy();      // Get initial energy
    fprintf(stderr, "Initial energy: %.10f\n", energy_ini);
    iterations = 0;

    //std::cout << "Crossing time: " << get_crossing_time() << std::endl;
    //std::cout << "Relaxation time: " << get_relaxation_time() << std::endl;
    //std::cout << "Kinetic energy: " << ekin << std::endl;
    //std::cout << "Potential energy: " << epot << std::endl;

    while (ITIME < int_time && iterations < ITE_MAX)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old();
        predicted_pos_vel(ITIME);
        update_acc_jerk(total);
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;

        if(iterations % (int)OUT == 0)
        {
            fprintf(stderr, "ITIME: %.10f ITE: %d\n", ITIME, iterations);
        }
        //get_energy_log(OUT, ITIME);
    }
    energy_end = energy();
    fprintf(stderr, "Final energy: %.10f\n", energy_end);
    fprintf(stderr, "DeltaE: %.10e\n", (energy_end - energy_ini)/energy_ini);
}
