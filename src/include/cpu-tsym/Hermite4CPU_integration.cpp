#include "Hermite4CPU.hpp"

void Hermite4CPU::integration()
{

    ns->gtime.integration_ini = omp_get_wtime();

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;

    // Setting maximum number of threads for OpenMP sections
    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    // Initial energy calculation
    ns->en.ini = nu->get_energy();
    ns->en.tmp = ns->en.ini;

    // Getting system information:
    // * Crossing Time and Half-mass Relaxation Time
    // * Close Encounter Radius and Timestep
    nu->nbody_attributes();

    update_neighbour_radius();

    init_acc_jrk(ns->h_p, ns->h_f, ns->h_r_sphere);
    init_dt(ATIME, ETA_S);

    logger->print_info();
    logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, ns->en.ini);

    if (ns->ops.print_all)
    {
        logger->print_all(ITIME);
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


        // If we have MultipleSystems already created
        // we procede with the time-symmetric integration
        if ((int)ms.size() > 0)
        {
          multiple_systems_integration(ms, ITIME);
          //printf("ITIME: %.15e | %.15e \n", ITIME, ns->integration_time);
        }

        // Considerar particulas virtuales y con masa cero
        predicted_pos_vel(ITIME, ns->h_t, ns->h_r, ns->h_v, ns->h_f, ns->h_p);
        update_acc_jrk(nact, ns->h_move, ns->h_r_sphere, ns->h_p, ns->h_f);

        // TODO: Check for encounters between single stars or binaries
        //print_nb(ITIME, nb_list, ns->h_f, ns->n, ns->h_p, ns->h_r_sphere);
        new_binaries = false;
        new_binaries = get_close_encounters(ITIME, nb_list, ns->h_f, ns->n, ns->h_p,
                                            ns->h_r_sphere, pairs, nact);


        correction_pos_vel(ITIME, nact, ns->h_dt, ns->h_t, ns->h_move,
                           ns->h_p, ns->h_f, ns->h_old, ns->h_a2, ns->h_a3,
                           ns->h_r, ns->h_v);

        // TODO: If an encounter in previous step
        // binary construct (replacing, virtual, etc)
        // passing the index of the system that it's created.
        // Consider the case that a new particle enters into the
        // system
        if(new_binaries)
        {
            for (int b = 0; b < (int)pairs.size(); b++)
            {
                // New MultipleSystem
                std::cout << "Creating new multiple system "
                          << pairs[b].id_a << " , " << pairs[b].id_b
                          << std::endl;

                MultipleSystem new_ms(ns, nu);

                int id_a = pairs[b].id_a;
                int id_b = pairs[b].id_b;

                // Adding the binary ids
                new_ms.add_particle(id_a);
                new_ms.add_particle(id_b);

                // Initialization of the binary
                new_ms.init_timestep();

                // ghost particle which will be store in the first member
                // of the new binary.
                create_ghost_particle(new_ms);

                // The second member of the binary will remain in the system
                // but its mass will be `0`, so in this way we avoid removing
                // this particle and moving all the system, which is computationally
                // expensive.
                // This particle will not affect the evolution of the system
                // since the force he will contribute is zero.

                // Adding the new binary to the vector
                ms.push_back(new_ms);
                pairs.erase(pairs.begin());
            }
        }

        // TODO: Termination
        // Termination of a simple binary occurs when the distance
        // between its members becomes greater than R cl .
        //  Hard binaries are not terminated, unless another particle
        // becomes a member and interacts strongly with their members.
        // Collisions ?

        if ((int)ms.size() > 0)
        {
            for (int i = 0; i < (int)ms.size(); i++)
            {
                double rx = ms[i].parts[1].r.x - ms[i].parts[0].r.x;
                double ry = ms[i].parts[1].r.y - ms[i].parts[0].r.y;
                double rz = ms[i].parts[1].r.z - ms[i].parts[0].r.z;

                double r = sqrt(rx * rx + ry * ry + rz * rz);
                if ( r > ns->r_cl)
                {
                    std::cout << "Termination!" << std::endl;
                    getchar();
                }
            }

        }


        // Update the amount of interactions counter
        interactions += nact * ns->n;

        // Find the next integration time
        next_integration_time(ATIME);

        if(std::ceil(ITIME) == ITIME)
        {
            assert(nact == ns->n);
            logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy());
            if (ns->ops.print_all)
            {
                logger->print_all(ITIME);
            }
            if (ns->ops.print_lagrange)
            {
                nu->lagrange_radii();
                logger->print_lagrange_radii(ITIME, nu->layers_radii);
            }
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        ns->iterations++;
    }
    ns->gtime.integration_end =  omp_get_wtime() - ns->gtime.integration_ini;
}
