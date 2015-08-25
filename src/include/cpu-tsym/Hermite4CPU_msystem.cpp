#include "Hermite4CPU.hpp"

void Hermite4CPU::multiple_systems_integration(std::vector<MultipleSystem> &ms, double ITIME, int **nb_list)
{
    // P(EC)^3
    for (int i = 0; i < (int)ms.size(); i++)
    {

        // ITIME = Final time
        // CTIME = Current time

        // The time of all the members will be the same
        double CTIME = ms[i].parts[0].t + ms[i].parts[0].dt;

        //ms[i].get_orbital_elements();
        //getchar();
        //ms[i].ini_e = ms[i].get_energy();

        long long int iterations = 0;
        while (CTIME < ITIME)
        {
            ////ms[i].get_orbital_elements();
            ms[i].prediction(CTIME);
            ms[i].save_old();

            // (EC)^3
            for (int k = 0; k < 3; k++)
            {
                ms[i].evaluation(nb_list[ms[i].parts[0].id]);
                ms[i].perturbers(CTIME, k);
                ms[i].correction(CTIME, true);
            }

            ms[i].update_timestep(CTIME);
            ms[i].next_itime(CTIME);
            iterations++;
            //printf("33 1 %.5e %.5e %.5e %.5e %.5e %.5e\n", ms[i].parts[0].r.x, ms[i].parts[0].r.y, ms[i].parts[0].r.z, ms[i].parts[0].v.x, ms[i].parts[0].v.y, ms[i].parts[0].v.z);
            //printf("33 2 %.5e %.5e %.5e %.5e %.5e %.5e\n", ms[i].parts[1].r.x, ms[i].parts[1].r.y, ms[i].parts[1].r.z, ms[i].parts[1].v.x, ms[i].parts[1].v.y, ms[i].parts[1].v.z);
        }
        //double end_e = ms[i].get_energy();
        //printf("End binary evolution: DE = %.15e E = %.15e | Ite: %lld\n", (end_e - ms[i].ini_e)/ms[i].ini_e, end_e, iterations);
    }
}

// Update the neighbour radius based on the close encounter radius
void Hermite4CPU::update_neighbour_radius()
{
    #pragma omp parallel for
    for (int i = 0; i < ns->n; i++)
    {
        // MYRIAD
        double mass = ns->h_r[i].w + ns->m_g;
        ns->h_r_sphere[i] = 50 * sqrt((ns->n * mass) * 0.5) * ns->r_cl;
        //ns->h_r_sphere[i] = 50 * sqrt((ns->n * mass) * 0.5) * ns->r_cl;
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
    //getchar();
}

bool Hermite4CPU::get_close_encounters(double itime, int **nb_list, Forces *f,
                                       int n, Predictor *p, double *r_sphere,
                                       std::vector<binary_id> &pairs, int nact)
{
    for (int k = 0; k < nact; k++)
    {
        int i = ns->h_move[k];
        // Amount of neighbours of the i-particle
        int nb = f[i].nb;
        //TODO: Check if i-particle is a ghost particle
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
                            //if (kin - pot < -0.01)
                            if (kin - pot < 0)
                            {
                                printf("Adding a pair %d %d | Binding energy %.5e\n", i, k, kin-pot);
                                //printf("Binary: r: %.5e r_crit: %.5e\n",r, r_crit);
                                //printf("Binary: k: %.5e k_sys: %.5e\n", kin, ns->en.kinetic);
                                //printf("Binary: dt: %.5e dt: %.5e | dt: %.5e\n", ns->h_dt[i], ns->h_dt[k], ns->dt_cl);

                                // This particles will dissapear in the next
                                // step (inside the integration loop), that is why
                                // they are now ghost.
                                ghosts[i] = 1;
                                ghosts[k] = 1;

                                // This `pairs` variable will contain the binary
                                // candidate.
                                binary_id bin = {i, k};
                                pairs.push_back(bin);

                                // Assigning minimum timestep to all the neighbors
                                //for (int z = 0; z < f[i].nb; z++)
                                //{
                                //    int n_id = nb_list[i][z];
                                //    if (n_id != k && n_id != i)
                                //    {
                                //        printf(">>dt min to %d\n", n_id);
                                //        ns->h_dt[n_id] = D_TIME_MIN;
                                //    }
                                //}
                            }
                        }
                    }
                    else
                    {
                        // TODO
                        // Second member is a ghost particle
                        // We can check the mass to see which one it's
                    }
                }
            }
        }
        else
        {
            // TODO
            // First member is a ghost particle
        }
    }

    if ((int)pairs.size() > 0)
        return true;

    return false;
}

SParticle Hermite4CPU::create_ghost_particle(MultipleSystem ms)
{
    // Getting center of mass of the new multiple system

    SParticle sp = ms.get_center_of_mass(ms.parts[0], ms.parts[1]);

    //printf("CoM %.15e %.15e %.15e\n", sp.r.x, sp.r.y, sp.r.z);

    // Replacing first member by a ghost particle
    int id = ms.parts[0].id;

    //std::cout << "Creating ghost particle between "
    //          << ms.parts[0].id << " and " << ms.parts[1].id
    //          << " (using id " << ms.parts[0].id << ")" << std::endl;


    ns->h_r[id]  = sp.r; // The mass is in the .w component
    ns->h_v[id]  = sp.v;
    ns->h_f[id]  = sp.f;
    ns->h_old[id]  = sp.old;

    ns->h_dt[id] = D_TIME_MIN;
//    ns->h_spehere[id]

    // Setting a zero mass to the second member.
    ns->h_r[ms.parts[1].id].w = 0.0;

    // Maybe avoid having only an empty particle active.
    ns->h_dt[ms.parts[1].id] = D_TIME_MAX;

    // TODO: Maybe add a blacklist to the method "find_particles_to_move"
    // using all the particles that have no mass.

    return sp;
}

