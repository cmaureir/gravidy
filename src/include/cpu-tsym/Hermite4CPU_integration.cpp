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
    double ms_energy = 0.0;

    bool dummy = false;

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
          multiple_systems_integration(ms, ITIME, nb_list);
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

                // ghost particle which will be store in the first member
                // of the new binary.
                SParticle sp = create_ghost_particle(new_ms);

                new_ms.adjust_particles(sp);

                // Initialization of the binary
                new_ms.evaluation(NULL);
                new_ms.init_timestep();


                dummy = true;
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

        // Update the amount of interactions counter
        interactions += nact * ns->n;

        // Find the next integration time
        next_integration_time(ATIME);

        if(std::ceil(ITIME) == ITIME || new_binaries || dummy)
        {
            //assert(nact == ns->n);

            // Check for MultipleSystems and get the energy.
            for (int i = 0; i < (int)ms.size(); i++)
            {
                ms_energy += ms[i].get_energy();
            }
            printf("MS Energy: %.15e\n", ms_energy);

            //logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy(ms_energy));
            logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy(0));
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

        // TODO: Termination, dectect hard binaries
        // Termination of a simple binary occurs when the distance
        // between its members becomes greater than R cl .
        // Hard binaries are not terminated, unless another particle
        // becomes a member and interacts strongly with their members.

        if ((int)ms.size() > 0)
        {
            for (int i = 0; i < (int)ms.size(); i++)
            {
                MParticle part0 = ms[i].parts[0];
                MParticle part1 = ms[i].parts[1];

                double rx = part1.r.x - part0.r.x;
                double ry = part1.r.y - part0.r.y;
                double rz = part1.r.z - part0.r.z;

                double r = sqrt(rx * rx + ry * ry + rz * rz);
                if ( r > ns->r_cl)
                {
                    std::cout << "Termination!" << std::endl;
                    int id0 = part0.id;
                    int id1 = part1.id;

                    //SParticle sp = ms[i].get_center_of_mass(part0, part1);

                    // Part1
                    // Part 1 = CoM particle + Current Part 1 position/velocity
                    ns->h_r[id1].x  = ns->h_r[id0].x + part1.r.x;
                    ns->h_r[id1].y  = ns->h_r[id0].y + part1.r.y;
                    ns->h_r[id1].z  = ns->h_r[id0].z + part1.r.z;

                    ns->h_v[id1].x  = ns->h_r[id0].x + part1.v.x;
                    ns->h_v[id1].y  = ns->h_r[id0].y + part1.v.y;
                    ns->h_v[id1].z  = ns->h_r[id0].z + part1.v.z;
                    ns->h_r[id1].w = part1.r.w;;

                    ns->h_f[id1].a[0]  = ns->h_f[id0].a[0] + part0.f.a[0];
                    ns->h_f[id1].a[1]  = ns->h_f[id0].a[1] + part0.f.a[1];
                    ns->h_f[id1].a[2]  = ns->h_f[id0].a[2] + part0.f.a[2];

                    ns->h_f[id1].a1[0] = ns->h_f[id0].a1[0] + part0.f.a1[0];
                    ns->h_f[id1].a1[1] = ns->h_f[id0].a1[1] + part0.f.a1[1];
                    ns->h_f[id1].a1[2] = ns->h_f[id0].a1[2] + part0.f.a1[2];

                    //ns->h_a2[id1] = ns->h_a2[id0] + part1.a2;
                    //ns->h_a3[id1] = ns->h_a3[id0] + part1.a3;
                    //ns->h_t[id1]  = ns->h_t[id0]  + part1.t;
                    ns->h_dt[id1] = D_TIME_MIN;

                    // Part0
                    ns->h_r[id0].x  += part0.r.x;
                    ns->h_r[id0].y  += part0.r.y;
                    ns->h_r[id0].z  += part0.r.z;
                    ns->h_r[id0].w = part0.r.w;;

                    ns->h_v[id0].x  += part0.v.x;
                    ns->h_v[id0].y  += part0.v.y;
                    ns->h_v[id0].z  += part0.v.z;

                    ns->h_f[id0].a[0]  += part0.f.a[0];
                    ns->h_f[id0].a[1]  += part0.f.a[1];
                    ns->h_f[id0].a[2]  += part0.f.a[2];

                    ns->h_f[id0].a1[0]  += part0.f.a1[0];
                    ns->h_f[id0].a1[1]  += part0.f.a1[1];
                    ns->h_f[id0].a1[2]  += part0.f.a1[2];

                    //ns->h_a2[id0].x += part0.a2.x;
                    //ns->h_a2[id0].y += part0.a2.y;
                    //ns->h_a2[id0].z += part0.a2.z;

                    //ns->h_a3[id0].x += part0.a3.x;
                    //ns->h_a3[id0].y += part0.a3.y;
                    //ns->h_a3[id0].z += part0.a3.z;

                    //ns->h_t[id0]  += part0.t;
                    ns->h_dt[id0] += D_TIME_MIN;


printf("New P1 %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", ns->h_r[id0].w, ns->h_r[id0].x, ns->h_r[id0].y, ns->h_r[id0].z, ns->h_v[id0].x, ns->h_v[id0].y, ns->h_v[id0].z);
printf("New P2 %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", ns->h_r[id1].w, ns->h_r[id1].x, ns->h_r[id1].y, ns->h_r[id1].z, ns->h_v[id1].x, ns->h_v[id1].y, ns->h_v[id1].z);

            ms_energy = ms[i].get_energy();
            logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy(ms_energy));

                    ms.erase(ms.begin()+i, ms.begin()+i+1);

                }
            }
        }
        ms_energy = 0.0;

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        ns->iterations++;
    }
    ns->gtime.integration_end =  omp_get_wtime() - ns->gtime.integration_ini;
}
