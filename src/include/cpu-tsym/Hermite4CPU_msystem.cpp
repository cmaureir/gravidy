#include "Hermite4CPU.hpp"

void Hermite4CPU::msystem_integration(std::vector<MultipleSystem> &ms, double ITIME, int **nb_list)
{
    // P(EC)^3
    for (int i = 0; i < (int)ms.size(); i++)
    {

        // ITIME = Final time
        // CTIME = Current time

        // The time of all the members will be the same
        double CTIME = ms[i].parts[0].t + ms[i].parts[0].dt;

        ms[i].get_orbital_elements(true);

        ms[i].ini_e = ms[i].get_energy();
        ms[i].evaluation(NULL);
        ms[i].perturbers(ms[i].parts[0].t, -1);

        long long int iterations = 0;
        while (CTIME < ITIME)
        {
            // (P) Prediction step
            ms[i].prediction(CTIME);

            ms[i].save_old();

            // (EC)^3
            for (int k = 0; k < 3; k++)
            {
                // (E) Evaluation of the forces step
                ms[i].evaluation(nb_list[ms[i].parts[0].id]);

                // Perturbers effect
                ms[i].perturbers(CTIME, k);

                // (C) Correction step
                ms[i].correction(CTIME, true);
            }

            ms[i].adjust_to_center_of_mass();
            ms[i].update_timestep(CTIME);
            ms[i].next_itime(CTIME);
            iterations++;
        }

        ms[i].get_orbital_elements(false);
        double end_e = ms[i].get_energy();
        printf("End binary evolution: DE = %.15e E = %.15e | ecc = %.5e a = %.5e | Ite: %lld\n",
                (end_e - ms[i].ini_e)/ms[i].ini_e,
                end_e,
                ms[i].ecc_end,
                ms[i].semimajor_end,
                iterations);
    }
}

void Hermite4CPU::msystem_termination(std::vector<MultipleSystem> &ms)
{

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

            double vx = part1.v.x - part0.v.x;
            double vy = part1.v.y - part0.v.y;
            double vz = part1.v.z - part0.v.z;

            //if ( r > ns->r_cl)
            double r_crit = sqrt(0.5 * ns->n * (part0.r.w + part1.r.w)) * ns->r_cl;

            if (r > r_crit && (rx*vx + ry*vy + rz*vz) > 0)
            {
                // The first particle ID represent the CoM (virtual)
                // particle in the system, so we need to treat it different.
                int id0 = part0.id;
                int id1 = part1.id;
                ms[i].get_orbital_elements(false);
                printf("05 Termination (%d, %d) EF: %.6e | a = %.4e | e = %.4e | pert = %d\n",
                    id0, id1, ms[i].get_energy(), ms[i].semimajor_end,
                    ms[i].ecc_end, ms[i].num_pert);

                // Part1
                // Part 1 = CoM particle + Current Part 1 position/velocity
                ns->h_r[id1] = ns->h_r[id0] + part1.r;
                ns->h_v[id1] = ns->h_v[id0] + part1.v;
                ns->h_r[id1].w = part1.r.w;;
                ns->h_f[id1] = ns->h_f[id0] + part1.f;
                ns->h_old[id1] = ns->h_f[id1];
                ns->h_t[id1] = ns->h_t[id0];
                ns->h_dt[id1] = D_TIME_MIN;

                // Part0
                ns->h_r[id0] += part0.r;
                ns->h_v[id0] += part0.v;
                ns->h_r[id0].w = part0.r.w;;
                ns->h_f[id0] += part0.f;
                ns->h_old[id0] = ns->h_f[id0];
                ns->h_t[id0] = ns->h_t[id0];
                ns->h_dt[id0] = D_TIME_MIN;

                virtuals[id0] = 0;
                ghosts[id1] = 0;

                ms.erase(ms.begin()+i, ms.begin()+i+1);

            }
        }
    }
}

// Update the neighbour radius based on the close encounter radius
void Hermite4CPU::update_neighbour_radius()
{
    //#pragma omp parallel for
    for (int i = 0; i < ns->n; i++)
    {
        // MYRIAD
        double mass = ns->h_r[i].w + ns->m_g;
        ns->h_r_sphere[i] = 5 * sqrt((ns->n * mass) * 0.5) * ns->r_cl;
    }
}

// Print all the neighbours for every particle
void Hermite4CPU::print_nb(double itime, int **nb_list, Forces *f, int n,
                           Predictor *p, double *r_sphere)
{
    for (int i = 0; i < n; i++)
    {
        int nb = f[i].nb;
        if (nb > 0)
        {
            std::cout << "Part[" << std::setw(2) << i << "]: ";
            for (int j = 0; j < nb; j++)
            {
                int k = nb_list[i][j];
                std::cout << k << " ";
            }
            std::cout << std::endl;
        }
    }
}

bool Hermite4CPU::get_close_encounters(double itime, int **nb_list, Forces *f,
                                       int n, Predictor *p, double *r_sphere,
                                       std::vector<binary_id> &pairs, int nact)
{

    // Checking close encounters just for active particles
    for (int a = 0; a < nact; a++)
    {
        int i = ns->h_move[a];

        // Amount of neighbours of the i-particle
        int nb = f[i].nb;

        // We don't consider ghost particles (m=0) since they are not important
        // for the system
        if(!ghosts[i])
        {
            if (nb > 0)
            {
                for (int j = 0; j < nb; j++)
                {
                    // k, Id of the j-neighbour
                    int k = nb_list[i][j];
                    if (i == k) continue;

                    if (!ghosts[k])
                    {
                        // Load the Predictor elements once, not everytime.
                        Predictor pk = p[k];
                        Predictor pi = p[i];

                        double rx = pk.r[0] - pi.r[0];
                        double ry = pk.r[1] - pi.r[1];
                        double rz = pk.r[2] - pi.r[2];

                        double r = sqrt(rx * rx + ry * ry + rz * rz);
                        double r_crit = sqrt(0.5 * ns->n * (pk.m + pi.m)) * ns->r_cl;

                        if (r <= r_crit)
                        {
                            double vx = pk.v[0] - pi.v[0];
                            double vy = pk.v[1] - pi.v[1];
                            double vz = pk.v[2] - pi.v[2];

                            double v2 = vx*vx + vy*vy + vz*vz;

                            double kin = 0.5 * v2;
                            double pot = (pk.m + pi.m)/r;

                            // Binding energy
                            if (kin - pot < 0 && (rx*vx + ry*vy + rz*vz) < 0)
                            {

                                if (i < k)
                                {

                                    if (virtuals[i])
                                    {
                                        printf("i-particle %d is a virtual particle (companion %d)\n", i, k);
                                        getchar();
                                    }
                                    else if (virtuals[k])
                                    {
                                        printf("k-particle %d is a virtual particle (companion %d)\n", k, i);
                                        getchar();
                                    }

                                    // First particle will be virtual and the second
                                    // a ghost particle.
                                    virtuals[i] = 1;
                                    ghosts[k] = 1;
                                    // New binary's ID
                                    binary_id bin = {i, k};
                                    pairs.push_back(bin);
                                    printf("New close encounter %d and %d\n", i, k);
                                }
                                else
                                {
                                    printf("This should never happen...");
                                    getchar();
                                }
                            }
                        }
                    } // if(!ghosts[k])
                }
            } // if (nb > 0)
        } // if (!ghosts[k])
    }

    if ((int)pairs.size() > 0)
        return true;

    return false;
}

SParticle Hermite4CPU::create_virtual_particle(MultipleSystem ms)
{
    // Getting center of mass of the new multiple system

    SParticle sp = ms.get_center_of_mass(ms.parts[0], ms.parts[1]);

    //printf("CoM %.15e %.15e %.15e\n", sp.r.x, sp.r.y, sp.r.z);

    // Replacing first member by a virtual particle
    int id = ms.parts[0].id;

    ns->h_r[id]  = sp.r; // The mass is in the .w component
    ns->h_v[id]  = sp.v;
    ns->h_f[id]  = sp.f;
    ns->h_old[id]  = sp.old;

    ns->h_dt[id] = D_TIME_MIN;
    //ns->h_spehere[id]

    // Setting a zero mass to the second member.
    ns->h_r[ms.parts[1].id].w = 0.0;

    // Maybe avoid having only an empty particle active.
    ns->h_dt[ms.parts[1].id] = D_TIME_MAX;

    // TODO: Maybe add a blacklist to the method "find_particles_to_move"
    // using all the particles that have no mass.

    return sp;
}


void Hermite4CPU::unpack_and_get_energy(std::vector<MultipleSystem> ms, double ITIME, long long int interactions, int nsteps)
{

    double4 tmp_r0[100];
    double4 tmp_r1[100];
    double4 tmp_v0[100];
    double4 tmp_v1[100];

    //logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy());
    for (int i = 0; i < (int)ms.size(); i++)
    {
        MParticle part0 = ms[i].parts[0];
        MParticle part1 = ms[i].parts[1];

        int id0 = part0.id;
        int id1 = part1.id;

        //printf("[%d] p1 % .9e | % .9e % .9e % .9e (% .9e % .9e % .9e)\n", i,
        //    ns->h_r[id0].w, ns->h_r[id0].x, ns->h_r[id0].y, ns->h_r[id0].z,
        //                    ns->h_v[id0].x, ns->h_v[id0].y, ns->h_v[id0].z);

        //printf("[%d] p2 % .9e | % .9e % .9e % .9e (% .9e % .9e % .9e)\n", i,
        //    ns->h_r[id1].w, ns->h_r[id1].x, ns->h_r[id1].y, ns->h_r[id1].z,
        //                    ns->h_v[id1].x, ns->h_v[id1].y, ns->h_v[id1].z);

        //printf("BB  p1 % .9e | % .9e % .9e % .9e (% .9e % .9e % .9e)\n",
        //    part0.r.w, part0.r.x, part0.r.y, part0.r.z,
        //               part0.v.x, part0.v.y, part0.v.z);
        //printf("BB  p2 % .9e | % .9e % .9e % .9e (% .9e % .9e % .9e)\n",
        //    part1.r.w, part1.r.x, part1.r.y, part1.r.z,
        //               part1.v.x, part1.v.y, part1.v.z);

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

        //printf("Moving\n");
        //printf("[%d] p1 % .9e | % .9e % .9e % .9e (% .9e % .9e % .9e)\n", i, ns->h_r[id0].w, ns->h_r[id0].x, ns->h_r[id0].y, ns->h_r[id0].z, ns->h_v[id0].x, ns->h_v[id0].y, ns->h_v[id0].z);
        //printf("[%d] p2 % .9e | % .9e % .9e % .9e (% .9e % .9e % .9e)\n", i, ns->h_r[id1].w, ns->h_r[id1].x, ns->h_r[id1].y, ns->h_r[id1].z, ns->h_v[id1].x, ns->h_v[id1].y, ns->h_v[id1].z);

    }

    logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy());


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
    //logger->print_energy_log(ITIME, ns->iterations, interactions, nsteps, nu->get_energy());
}
