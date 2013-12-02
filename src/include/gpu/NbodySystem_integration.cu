#include "../NbodySystem.hpp"

void NbodySystem::integration(Hermite4GPU h, Logger log, NbodyUtils nu)
{

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;


    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    h.set_pointers(d_p, d_i, h_i, d_fout, d_fout_tmp, h_fout_tmp, d_f, d_move);

    h.init_acc_jrk(h_p, h_f);
    h.init_dt(ATIME, h_dt, h_t, h_f);

    en.ini = get_energy_gpu();   // Initial calculation of the energy of the system
    en.tmp = en.ini;

    nu.set_energy(en.tmp);

    hmr_time = nu.get_half_mass_relaxation_time();
    cr_time  = nu.get_crossing_time();

    log.print_info(n, e2, eta, integration_time, hmr_time, cr_time);
    log.print_energy_log(ITIME, iterations, interactions, nsteps, gtime, en, en.ini);

    if (ops.print_all)
    {
        log.print_all(ITIME, n, h_r, h_v, h_f, h_dt);
    }
    if (ops.print_lagrange)
    {
        nu.lagrange_radii();
        log.print_lagrange_radii(ITIME, nu.layers_radii);
    }

    while (ITIME < integration_time)
    {
        ITIME = ATIME;

        nact = h.find_particles_to_move(h_move, ITIME, h_dt, h_t);

        h.save_old_acc_jrk(nact, h_move, h_old, h_f);

        h.predicted_pos_vel(ITIME, h_p, h_r, h_v, h_f, h_t, gtime);

        h.update_acc_jrk(nact, h_move, h_p, h_f, gtime);

        h.correction_pos_vel(ITIME, nact, h_move, h_r, h_v, h_f, h_t, h_dt, h_p, h_old, h_a3, h_a2, gtime);

        // Update the amount of interactions counter
        interactions += nact * n;

        // Find the next integration time
        h.next_integration_time(ATIME, h_dt, h_t);


        if(std::ceil(ITIME) == ITIME)
        //if(nact == n)
        {
            assert(nact == n);
            log.print_energy_log(ITIME, iterations, interactions, nsteps, gtime, en, get_energy());
            if (ops.print_all)
            {
                log.print_all(ITIME, n, h_r, h_v, h_f, h_dt);
            }
            if (ops.print_lagrange)
            {
                nu.lagrange_radii();
                log.print_lagrange_radii(ITIME, nu.layers_radii);
            }
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        iterations++;
    }
}
