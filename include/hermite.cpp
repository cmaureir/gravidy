#include "hermite.hpp"


void integrate_gpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int total;
    iterations = 0;

    gpu_init_acc_jrk();
    init_dt(&ATIME);

    energy_ini = gpu_energy();      // Get initial energy

    while (ITIME < int_time)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old(total);
        gpu_predicted_pos_vel(ITIME);
        gpu_update_acc_jrk(total);
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;
        nsteps += total;

        if(std::ceil(ITIME) == ITIME)
        {
            get_energy_log(ITIME, iterations, nsteps);
        }

    }
}

void integrate_cpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int total = 0;
    int nsteps = 0;
    iterations = 0;

    init_acc_jrk();
    init_dt(&ATIME);
    energy_ini = energy();
    energy_tmp = 0.0;

//    Point p = get_center_of_density();
//    for (int i = 0; i < n; i++) {
//        double rx = h_r[i].x - p.x;
//        double ry = h_r[i].y - p.y;
//        double rz = h_r[i].z - p.z;
//        double d  = get_magnitude(rx, ry, rz);
//        printf("%d %.6f\n",i, d);
//    }

//    float t_rh = get_relaxation_time();
//    float t_cr = get_crossing_time();
//
//    std::cout << "T_rh : " << t_rh << std::endl;
//    std::cout << "T_cr : " << t_cr << std::endl;
//    std::cout << "T_cc : " << 17 * t_rh << std::endl;
//
//    printf("%d\n",n);
    while (ITIME < int_time)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old(total);
        predicted_pos_vel(ITIME);
        update_acc_jrk(total);
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;
        nsteps += total;

        if(std::ceil(ITIME) == ITIME)
        {
            get_energy_log(ITIME, iterations, nsteps);
        }

    }
}
