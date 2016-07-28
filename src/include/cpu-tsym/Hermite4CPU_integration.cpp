#include "Hermite4CPU.hpp"

void Hermite4CPU::integration()
{

    ns->gtime.integration_ini = omp_get_wtime();

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = ns->snapshot_time;     // Integration time
    unsigned int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;

    // Setting maximum number of threads for OpenMP sections
    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    update_neighbour_radius();

    init_acc_jrk(ns->h_p, ns->h_f, ns->h_r_sphere);
    init_dt(ATIME, ETA_S, ITIME);

    // Initial energy calculation
    ns->en.ini = nu->get_energy();
    ns->en.tmp = ns->en.ini;

    // Getting system information:
    nu->nbody_attributes();

    //update_neighbour_radius();


    logger->print_info();
    logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, ns->en.ini);

    int snap_number = ns->snapshot_number;
    logger->write_snapshot(snap_number, ITIME);
    snap_number++;

    if (ns->ops.print_all)
    {
        logger->print_all(ITIME, snap_number);
    }
    if (ns->ops.print_lagrange)
    {
        nu->lagrange_radii();
        logger->print_lagrange_radii(ITIME, nu->layers_radii);
    }


    // TODO: Move this variables from here
    bool new_binaries;
    std::vector<binary_id> pairs;
    std::vector<MultipleSystem> ms;

    while (ITIME < ns->integration_time)
    {
        // Current integration time
        ITIME = ATIME;

        nact = find_particles_to_move(ITIME);

        save_old_acc_jrk(nact);

        // MultipleSystems time-symmetric integration
        if ((int)ms.size())
        {
            msystem_integration(ms, ITIME, nb_list);
        }


        // (P) Prediction step
        predicted_pos_vel(ITIME, ns->h_t, ns->h_r, ns->h_v, ns->h_f, ns->h_p);

        // (E) Forces evaluation step
        update_acc_jrk(nact, ns->h_move, ns->h_r_sphere, ns->h_p, ns->h_f);

        //print_nb(ITIME, nb_list, ns->h_f, ns->n, ns->h_p, ns->h_r_sphere);

        new_binaries = false;
        new_binaries = get_close_encounters(ITIME, nb_list, ns->h_f, ns->n, ns->h_p,
                                            ns->h_r_sphere, pairs, nact);


        // (C) Correction step
        correction_pos_vel(ITIME, nact, ns->h_dt, ns->h_t, ns->h_move,
                           ns->h_p, ns->h_f, ns->h_old, ns->h_a2, ns->h_a3,
                           ns->h_r, ns->h_v);

        // Binary creation
        if(new_binaries)
        {
            for (int b = 0; b < (int)pairs.size(); b++)
            {
                MultipleSystem new_ms(ns, nu);

                int id_a = pairs[b].id_a;
                int id_b = pairs[b].id_b;

                // Adding the binary ids
                new_ms.add_particle(id_a);
                new_ms.add_particle(id_b);

                SParticle sp = create_virtual_particle(new_ms);
                new_ms.adjust_particles(sp);

                // Initialization of the binary
                new_ms.init_timestep();
                new_ms.ini_e = new_ms.get_energy();
                new_ms.get_orbital_elements(true);

                printf("04 New MS (%d, %d) | E0 = %.6e | a = %.4e | ecc = %.4e | pert = %d\n",
                    pairs[b].id_a, pairs[b].id_b,
                    new_ms.ini_e, new_ms.semimajor_ini, new_ms.ecc_ini, new_ms.num_pert);

                // Adding the new binary to the vector
                ms.push_back(new_ms);
                pairs.erase(pairs.begin());
            }
        }

        // Update the amount of interactions counter
        interactions += nact * ns->n;

        // Find the next integration time
        next_integration_time(ATIME);

        if(nact == ns->n && ITIME >= D_TIME_MAX)
        {
            //assert(nact == ns->n);

            unpack_and_get_energy(ms, ITIME, interactions, nsteps);

            if (ns->ops.print_all)
            {
                logger->print_all(ITIME, snap_number);
            }
            if (ns->ops.print_lagrange)
            {
                nu->lagrange_radii();
                logger->print_lagrange_radii(ITIME, nu->layers_radii);
            }

            update_neighbour_radius();
        }

        // Check for termination of multiple systems
        msystem_termination(ms);

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        ns->iterations++;

    }
    ns->gtime.integration_end =  omp_get_wtime() - ns->gtime.integration_ini;
}
