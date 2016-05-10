#include "Hermite4GPU.cuh"

void Hermite4GPU::integration()
{
    ns->gtime.integration_ini = omp_get_wtime();

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = ns->snapshot_time;     // Integration time
    unsigned int nact     = 0;       // Active particles
    unsigned int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;

    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    init_acc_jrk();
    init_dt(ATIME, ETA_S, ITIME);
    //logger->print_all(0);
    //getchar();

    ns->en.ini = get_energy_gpu();   // Initial calculation of the energy of the system
    ns->en.tmp = ns->en.ini;

    //nu->nbody_attributes();

    logger->print_info();
    logger->write_info();
    logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, ns->en.ini);

    int snap_number = ns->snapshot_number;
    logger->write_snapshot(snap_number, ITIME);
    snap_number++;

    if (ns->ops.print_all)
    {
        logger->print_all(ITIME);
    }
    if (ns->ops.print_lagrange)
    {
        nu->lagrange_radii();
        logger->print_lagrange_radii(ITIME, nu->layers_radii);
    }

    while (ITIME < ns->integration_time)
    {
        ITIME = ATIME;

        nact = find_particles_to_move(ITIME);

        save_old_acc_jrk(nact);

        predicted_pos_vel(ITIME);

        update_acc_jrk(nact);

        correction_pos_vel(ITIME, nact);

        // Update the amount of interactions counter
        interactions += nact * ns->n;

        // Find the next integration time
        next_integration_time(ATIME);


        if(nact == ns->n)
        {
            //assert(nact == ns->n);
            logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, get_energy_gpu());
            if (ns->ops.print_all)
            {
                logger->print_all(ITIME);
            }
            if (ns->ops.print_lagrange)
            {
                nu->lagrange_radii();
                logger->print_lagrange_radii(ITIME, nu->layers_radii);
            }
            logger->write_snapshot(snap_number, ITIME);
            snap_number++;
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        ns->iterations++;

    }
    ns->gtime.integration_end =  omp_get_wtime() - ns->gtime.integration_ini;
    logger->write_snapshot(snap_number, ITIME);
    //logger->add_info(std::string("SnapshotNumber:"), std::to_string(snap_number));
    logger->add_info(std::string("SnapshotNumber:"), std::string(SSTR(snap_number)));
}
