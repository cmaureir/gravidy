#include "hermite.hpp"
#define ITE_MAX (1e6)
#define OUT (1e3)


void integrate_gpu()
{

    float ATIME = 1.0e+10; // Actual integration time
    float ITIME = 0.0;     // Integration time
    int total;

    // Initial Acceleration and Jerk calculation on CPU
    init_time = omp_get_wtime();
    gpu_init_acc_jrk();
    init_time = omp_get_wtime() - init_time;

    // Initialazing block time steps
    gpu_init_dt(&ATIME);

    energy_ini = gpu_energy();      // Get initial energy
    fprintf(stderr, "Initial energy: %.10f\n", energy_ini);
    iterations = 0;

    while (ITIME < int_time && iterations < ITE_MAX)
    //while (ITIME < int_time)
    {
        // Setting actual integration time
        ITIME = ATIME;

        // Save old values of the Acceleration and Jerk
        save_old();
        total = gpu_find_particles_to_move(ITIME);
        gpu_predicted_pos_vel(ITIME);
        //print_predicted(n);
//        if (total < BSIZE * 2)
//            gpu_update_acc_jrk_single(total);
//        else
//            gpu_update_acc_jrk(total);
            gpu_update_acc_jrk(total);
//        gpu_update_2d(total);
//        print_accelerations(n);
        gpu_correction_pos_vel(ITIME, total);
    //    get_energy_log(OUT, ITIME);
        gpu_next_itime(&ATIME);
        iterations++;

    }
    //energy_end = gpu_energy();
    //fprintf(stderr, "Final energy: %.10f\n", energy_end);
    //fprintf(stderr, "DeltaE: %.10e\n", (energy_end - energy_ini)/energy_ini);
}

void
integrate_cpu()
{
    float ATIME = 1.0e+10; // Actual integration time
    float ITIME = 0.0;     // Integration time
    int total = 0;

    init_time = omp_get_wtime();
    init_acc_jrk();
    init_time = omp_get_wtime() - init_time;
    //printf("Init_acc_jrk %.10f\n", init_time);

    init_dt(&ATIME);
    energy_ini = energy();
    energy_tmp = energy_ini;
    //fprintf(stderr, "Initial energy: %.10f\n", energy_ini);
    iterations = 0;
    //std::cout << "Crossing time: " << get_crossing_time() << std::endl;
    //std::cout << "Relaxation time: " << get_relaxation_time() << std::endl;
    //std::cout << "Kinetic energy: " << ekin << std::endl;
    //std::cout << "Potential energy: " << epot << std::endl;
    float out_param = 0.1;

    //while (ITIME < int_time && iterations < ITE_MAX)
    int nsteps = 0;
    while (ITIME < int_time)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old();
        //predicted_pos_vel_kepler(ITIME, total);
        predicted_pos_vel(ITIME);
        //print_predicted(n);
        update_acc_jrk(total);
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;
        nsteps += total;

        get_energy_log(OUT, ITIME,nsteps, &out_param);

    }
    //energy_end = energy();
    //fprintf(stderr, "ETA_N: %.10f\n", ETA_N);
    //fprintf(stderr, "SQR ETA_N: %.10f\n", sqrt(ETA_N));
    //fprintf(stderr, "Iterations: %d\n", iterations);
    //fprintf(stderr, "Final energy: %.10f\n", energy_end);
    //fprintf(stderr, "DeltaE: %.10e\n", (energy_end - energy_ini)/energy_ini);
    //fprintf(stderr, "Init force time (s) : %.10e\n", init_time);
}
