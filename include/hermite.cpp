#include "hermite.hpp"
#define ITE_MAX (1e6)
#define OUT (1e2)

//void
//integrate_gpu()
//{
//
//    double ATIME = 1.0e+10; // Actual integration time
//    double ITIME = 0.0;     // Integration time
//    int total;
//    int cpu_simple = 0 , gpu_single=0, gpu_simple=0;
//
//    double total_e = 0.0f;
//    double ee = 0.0f;
//
//    // Initial Acceleration and Jerk calculation on CPU
//    gpu_init_acc_jerk();
//
//    // Initialazing block time steps
//    init_dt(&ATIME);
//
//    energy_ini = gpu_energy(0);      // Get initial energy
//    ee = energy_ini;
//
//    iterations = 0;
//    while (ITIME < int_time && iterations < ITE_MAX)
//    {
//        ITIME = ATIME;
//        total = find_particles_to_move(ITIME);
//        save_old();
//        predicted_pos_vel(ITIME);
//        if (total < BSIZE && n < 1000)
//        {
//           update_acc_jerk(total);
//           cpu_simple++;
//        }
//        else if(total < BSIZE * 4)
//        {
//          gpu_update_acc_jerk_single(total);
//          gpu_single++;
//        }
//        else
//        {
//            gpu_update_acc_jerk_tile(total);
//            gpu_simple++;
//        }
//        correction_pos_vel(ITIME, total);
//        next_itime(&ATIME);
//        iterations++;
//
//        energy_end = gpu_energy(1);
//        total_e += abs(energy_end - energy_ini)/energy_ini;
//
//        if(iterations % OUT == 0)
//        {
//            printf("%6d %.10f %4d %.10e\n", iterations, ITIME, total, total_e/OUT);
//            total_e = 0.0f;
//            energy_ini = energy_end;
//        }
//    }
//    energy_end = gpu_energy(1);      // Get final energy
//    std::cout << "n: " << n << std::endl;
//    std::cout << "CPU update: " << cpu_simple << std::endl;
//    std::cout << "GPU single: " << gpu_single << std::endl;
//    std::cout << "GPU simple: " << gpu_simple << std::endl;
//    std::cout << "DeltaE: " << (energy_end - ee)/ee << std::endl;
//}

void
integrate_cpu()
{
    float ATIME = 1.0e+10; // Actual integration time
    float ITIME = 0.0;     // Integration time
    int total;

    init_acc_jerk();

    init_dt(&ATIME);
    energy_ini = initial_energy();      // Get initial energy
    float ee = energy_ini;
    energy_total = 0.0f;
    // Using half-mass radius
    //t_rh = get_relaxation_time();
    // Using virial radius
    //t_cr = get_crossing_time();
    iterations = 0;
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
    std::cout << "DeltaE: " << (energy_end - ee)/ee << std::endl;
}
