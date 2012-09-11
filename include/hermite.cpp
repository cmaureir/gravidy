#include "hermite.hpp"
#define ITE_MAX (500)
#define OUT (1)

void
integrate_gpu()
{

    float ATIME = 1.0e+10; // Actual integration time
    float ITIME = 0.0;     // Integration time
    int total;

    energy_total = 0.0f;
    float ee = 0.0f;

    // Initial Acceleration and Jerk calculation on CPU
    gpu_init_acc_jerk();

    // Initialazing block time steps
    init_dt(&ATIME);

    energy_ini = gpu_energy(0);      // Get initial energy
    ee = energy_ini;
    fprintf(stderr, "Initial energy: %.5f\n", ee);
    iterations = 0;

    print_all(n);
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
            gpu_update_acc_jerk_simple(total);
        }
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;

        energy_end = gpu_energy(1);      // Get final energy
        energy_total += abs(energy_end - energy_ini)/energy_ini;
        get_energy_log(OUT, ITIME);

    }
    fprintf(stderr, "DeltaE: %.5f\n", (energy_end - ee)/ee);
}

void
integrate_cpu()
{
    float ATIME = 1.0e+10; // Actual integration time
    float ITIME = 0.0;     // Integration time
    int total;

    energy_total = 0.0f;
    float ee = 0.0f;

    init_acc_jerk();

    init_dt(&ATIME);
    energy_ini = initial_energy();      // Get initial energy
    ee = energy_ini;
    fprintf(stderr, "Initial energy: %.5f\n", ee);
    iterations = 0;

    print_all(n);
    while (ITIME < int_time && iterations <  ITE_MAX)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old();
        predicted_pos_vel(ITIME);
        update_acc_jerk(total);
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;

        energy_end = energy();      // Get final energy
        energy_total += abs(energy_end - energy_ini)/energy_ini;
        get_energy_log(OUT, ITIME);
    }
    fprintf(stderr, "DeltaE: %.5f\n", (energy_end - ee)/ee);
}
