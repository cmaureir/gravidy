#include "hermite.hpp"


void integrate_gpu()
{
    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int total;
    int nsteps = 0;
    iterations = 0;

    CUDA_SAFE_CALL(cudaMemcpy(d_r,  h_r,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v,  h_v,  d4_size,cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_m,  h_m,  f1_size,cudaMemcpyHostToDevice));
    gpu_init_acc_jrk();
    CUDA_SAFE_CALL(cudaMemcpy(h_a,  d_a,  d4_size,cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_j,  d_j,  d4_size,cudaMemcpyDeviceToHost));
    init_dt(&ATIME);
    energy_ini = gpu_energy();      // Get initial energy
    energy_tmp = 0.0;
    printf("Energy_ini: %.15e\n", energy_ini);
//    print_all(n,0);


    CUDA_SAFE_CALL(cudaMemcpy(d_t,  h_t,  d1_size,cudaMemcpyHostToDevice));
    while (ITIME < int_time)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old(total);


        gpu_predicted_pos_vel(ITIME);
        //gpu_update_acc_jrk(total);
        gpu_update_acc_jrk_simple(total);
        //update_acc_jrk(total);
        //gpu_update_2d(total);
        gpu_correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;
        nsteps += total;


        if(std::ceil(ITIME) == ITIME)
        {
            gpu_get_energy_log(ITIME, iterations, nsteps);
        }
//        print_all(n, 0);
//        printf("%.10f %4d\n", ITIME, total);
//        if(iterations == 10005)
//            break;
    }
    energy_end = gpu_energy();
    printf("Energy_end: %.15e\n", energy_end);
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
    printf("Energy_ini: %.15e\n", energy_ini);
//    print_all(n,0);

    while (ITIME < int_time)
    {
        ITIME = ATIME;
        total = find_particles_to_move(ITIME);
        save_old(total);
        predicted_pos_vel(ITIME);
//        predicted_pos_vel_kepler(ITIME, total);
        update_acc_jrk(total);
        correction_pos_vel(ITIME, total);
        next_itime(&ATIME);
        iterations++;
        nsteps += total;

        if(std::ceil(ITIME) == ITIME)
        {
            get_energy_log(ITIME, iterations, nsteps);
        }
//        print_all(n, 0);
//        printf("%.10f %4d\n", ITIME, total);
//        if(iterations == 1000)
//            break;
    }
    energy_end = energy();
    printf("Energy_end: %.15e\n", energy_end);

}
