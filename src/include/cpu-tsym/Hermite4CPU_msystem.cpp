#include "Hermite4CPU.hpp"

void Hermite4CPU::multiple_systems_integration(std::vector<MultipleSystem> ms, double ITIME)
{

    std::cout << "Starting binary evolution..." << (int)ms.size() << std::endl;

    // P(EC)^3
    for (int i = 0; i < (int)ms.size(); i++)
    {
        // ITIME = Final time
        // CTIME = Current time
        double CTIME = 1e6;

        ms[i].evaluation();
        ms[i].init_timestep(CTIME);

        std::cout << "CTIME = " << CTIME << std::endl;

        double ini_e = ms[i].get_energy();
        double end_e;

        std::cout << "E0 = " << ini_e << std::endl;
        printf("00 %.15f %.15f %.15f %.15f %.15f\n",
                ms[i].parts[0].t,
                ms[i].parts[0].dt,
                ms[i].parts[0].r.x,
                ms[i].parts[0].r.y,
                ms[i].parts[0].r.z);

        printf("01 %.15f %.15f %.15f %.15f %.15f\n",
                ms[i].parts[1].t,
                ms[i].parts[1].dt,
                ms[i].parts[1].r.x,
                ms[i].parts[1].r.y,
                ms[i].parts[1].r.z);
        getchar();

        long long int iterations = 0;
        while (CTIME < ITIME)
        {
            //ms[i].get_orbital_elements();
            ms[i].prediction(CTIME);
            ms[i].save_old();

            // (EC)^3
            for (int k = 0; k < 1; k++)
            {
                ms[i].evaluation();
                ms[i].correction(CTIME);
            }

            ms[i].update_information(CTIME);
            //if(std::ceil(CTIME) == CTIME)
            //{
            //    // orbital elements, energy
            //}
            // Update current time t_c = t_c + dt_b (binary)
            // if t_c is not t_f, we repeat the process.

            //ms[i].update_timestep(CTIME);
            ms[i].next_itime(CTIME);
            //ms[i].print();
            iterations++;
        }
            end_e = ms[i].get_energy();
            printf("Adios: DE = %.15e | Ite: %lld\n", (end_e - ini_e)/ini_e, iterations);

    }
}

// Update the neighbour radius based on the close encounter radius
void Hermite4CPU::update_neighbour_radius()
{
    #pragma omp parallel for
    for (int i = 0; i < ns->n; i++)
    {
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
                            if (kin - pot < 0)
                            {

                                std::cout << "Adding a pair "
                                          << i << " " << k
                                          << std::endl;
                                // This particles will dissapear in the next
                                // step (inside the integration loop), that is why
                                // they are now ghost.
                                ghosts[i] = 1;
                                ghosts[k] = 1;

                                // This `pairs` variable will contain the binary
                                // candidate.
                                binary_id bin = {i, k};
                                pairs.push_back(bin);
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

void Hermite4CPU::create_ghost_particle(MultipleSystem ms)
{
    // Getting center of mass of the new multiple system
    SParticle sp = ms.get_center_of_mass(ms.parts[0], ms.parts[1]);

    // Replacing first member by a ghost particle
    int id = ms.parts[0].id;

    std::cout << "Creating ghost particle between "
              << ms.parts[0].id << " and " << ms.parts[1].id
              << " (using id " << ms.parts[0].id << ")" << std::endl;


    ns->h_r[id]  = sp.r; // The mass is in the .w component
    ns->h_v[id]  = sp.v;
    ns->h_f[id]  = sp.f;

    ns->h_dt[id] = D_TIME_MIN;

    // Setting a zero mass to the second member.
    ns->h_r[ms.parts[1].id].w = 0.0;
    // TODO: Maybe add a blacklist to the method "find_particles_to_move"
    // using all the particles that have no mass.
}

