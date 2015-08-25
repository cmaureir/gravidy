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

    //double ekin_ini = nu->get_kinetic();
    //double epot_ini = nu->get_potential();

    // Getting system information:
    // * Crossing Time and Half-mass Relaxation Time
    // * Close Encounter Radius and Timestep
    nu->nbody_attributes();
    std::cout << "r_virial = " << ns->r_virial << std::endl;
    std::cout << "r_cl     = " << ns->r_cl     << std::endl;

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
    //double ms_energy = 0.0;

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
        new_binaries = false;
        new_binaries = get_close_encounters(ITIME, nb_list, ns->h_f, ns->n, ns->h_p,
                                            ns->h_r_sphere, pairs, nact);


        correction_pos_vel(ITIME, nact, ns->h_dt, ns->h_t, ns->h_move,
                           ns->h_p, ns->h_f, ns->h_old, ns->h_a2, ns->h_a3,
                           ns->h_r, ns->h_v);

        // Binary creation
        //new_binaries = false;
        if(new_binaries)
        {
            //print_nb(ITIME, nb_list, ns->h_f, ns->n, ns->h_p, ns->h_r_sphere);
            //getchar();
            for (int b = 0; b < (int)pairs.size(); b++)
            {
                MultipleSystem new_ms(ns, nu);

                int id_a = pairs[b].id_a;
                int id_b = pairs[b].id_b;

                // Adding the binary ids
                new_ms.add_particle(id_a);
                new_ms.add_particle(id_b);

                SParticle sp = create_ghost_particle(new_ms);
                new_ms.adjust_particles(sp);

                // Initialization of the binary
                new_ms.init_timestep();
                //new_ms.ini_e = new_ms.get_energy();

                printf("> New MS (%d, %d) | E0 = %.6e\n", pairs[b].id_a, pairs[b].id_b, new_ms.ini_e);

                // Adding the new binary to the vector
                ms.push_back(new_ms);
                pairs.erase(pairs.begin());
            }
        }

        // Update the amount of interactions counter
        interactions += nact * ns->n;

        // Find the next integration time
        next_integration_time(ATIME);

        if(nact == ns->n)
        {
            //assert(nact == ns->n);

            double4 tmp_r0[100];
            double4 tmp_r1[100];
            double4 tmp_v0[100];
            double4 tmp_v1[100];

            //logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy(ms_energy));
            for (int i = 0; i < (int)ms.size(); i++)
            {
                MParticle part0 = ms[i].parts[0];
                MParticle part1 = ms[i].parts[1];

                int id0 = part0.id;
                int id1 = part1.id;

                printf("[%d] p1 %.5e | %.5e %.5e %.5e (%.5e %.5e %.5e)\n", i, ns->h_r[id0].w, ns->h_r[id0].x, ns->h_r[id0].y, ns->h_r[id0].z, ns->h_v[id0].x, ns->h_v[id0].y, ns->h_v[id0].z);
                printf("[%d] p2 %.5e | %.5e %.5e %.5e (%.5e %.5e %.5e)\n", i, ns->h_r[id1].w, ns->h_r[id1].x, ns->h_r[id1].y, ns->h_r[id1].z, ns->h_v[id1].x, ns->h_v[id1].y, ns->h_v[id1].z);

                printf("BB p1 %.9e | %.9e %.9e %.9e (%.9e %.9e %.9e)\n", part0.r.w, part0.r.x, part0.r.y, part0.r.z, part0.v.x, part0.v.y, part0.v.z);
                printf("BB p2 %.9e | %.9e %.9e %.9e (%.9e %.9e %.9e)\n", part1.r.w, part1.r.x, part1.r.y, part1.r.z, part1.v.x, part1.v.y, part1.v.z);

                tmp_r0[i] = ns->h_r[id0];
                tmp_v0[i] = ns->h_v[id0];

                tmp_r1[i] = ns->h_r[id1];
                tmp_v1[i] = ns->h_v[id1];

                ns->h_r[id1] = ns->h_r[id0] + part1.r;
                ns->h_v[id1] = ns->h_v[id0] + part1.v;
                ns->h_r[id1].w = part1.r.w;;

                ns->h_r[id0] += part0.r;
                ns->h_v[id0] += part0.v;
                ns->h_r[id0].w = part0.r.w;;
                printf("Moving\n");
                printf("[%d] p1 %.5e | %.5e %.5e %.5e (%.5e %.5e %.5e)\n", i, ns->h_r[id0].w, ns->h_r[id0].x, ns->h_r[id0].y, ns->h_r[id0].z, ns->h_v[id0].x, ns->h_v[id0].y, ns->h_v[id0].z);
                printf("[%d] p2 %.5e | %.5e %.5e %.5e (%.5e %.5e %.5e)\n", i, ns->h_r[id1].w, ns->h_r[id1].x, ns->h_r[id1].y, ns->h_r[id1].z, ns->h_v[id1].x, ns->h_v[id1].y, ns->h_v[id1].z);

            }

            logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy(0));

            for (int i = 0; i < (int)ms.size(); i++)
            {
                MParticle part0 = ms[i].parts[0];
                MParticle part1 = ms[i].parts[1];
                int id0 = part0.id;
                int id1 = part1.id;

                ns->h_r[id0] = tmp_r0[i];
                ns->h_v[id0] = tmp_v0[i];

                ns->h_r[id1] = tmp_r1[i];
                ns->h_v[id1] = tmp_v1[i];

            }

            if (ns->ops.print_all)
            {
                logger->print_all(ITIME);
            }
            if (ns->ops.print_lagrange)
            {
                nu->lagrange_radii();
                logger->print_lagrange_radii(ITIME, nu->layers_radii);
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
                        // The first particle ID represent the CoM (ghost)
                        // particle in the system, so we need to treat it different.
                        int id0 = part0.id;
                        int id1 = part1.id;
                        printf("Termination! (%d, %d) EF: %.6e\n", id0, id1, ms[i].get_energy());


                        // Part1
                        // Part 1 = CoM particle + Current Part 1 position/velocity
                        ns->h_r[id1] = ns->h_r[id0] + part1.r;
                        ns->h_v[id1] = ns->h_v[id0] + part1.v;
                        ns->h_r[id1].w = part1.r.w;;
                        ns->h_f[id1] = ns->h_f[id0] + part1.f;
                        ns->h_old[id1] = ns->h_f[id1];
                        ns->h_t[id1] = ns->h_t[id0];
                        //ns->h_dt[id1] = ns->h_dt[id0];
                        ns->h_dt[id1] = D_TIME_MIN;

                        // Part0
                        ns->h_r[id0] += part0.r;
                        ns->h_v[id0] += part0.v;
                        ns->h_r[id0].w = part0.r.w;;
                        ns->h_f[id0] += part0.f;
                        ns->h_old[id0] = ns->h_f[id0];
                        ns->h_t[id0] = ns->h_t[id0];
                        //ns->h_dt[id0] = ns->h_dt[id0];
                        ns->h_dt[id0] = D_TIME_MIN;
                        //logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy(0));

                        ms.erase(ms.begin()+i, ms.begin()+i+1);

                    }
                }
            }

            update_neighbour_radius();
        }

        // Setting binary energy to zero
        //ms_energy = 0.0;

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        ns->iterations++;
    }
    ns->gtime.integration_end =  omp_get_wtime() - ns->gtime.integration_ini;
}
