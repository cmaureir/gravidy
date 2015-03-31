#include "Hermite4CPU.hpp"

Hermite4CPU::Hermite4CPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
                         : Hermite4(ns, logger, nu)
{

    nb_list  = new int*[ns->n];
    //nb_number = 2 * sqrt(ns->n);
    nb_number = NBMAX;

    ghosts = new int[ns->n];
    std::fill(ghosts, ghosts + ns->n, 0);

    #pragma omp parallel for
    for(int i = 0; i < ns->n; ++i)
        nb_list[i] = new int[nb_number];

}

Hermite4CPU::~Hermite4CPU()
{
    for(int i = 0; i < ns->n; ++i)
        delete nb_list[i];

    delete nb_list;
    delete ghosts;
}

void Hermite4CPU::force_calculation(Predictor pi, Predictor pj, Forces &fi,
                                     int i, int j, double hi)
{
    double rx = pj.r[0] - pi.r[0];
    double ry = pj.r[1] - pi.r[1];
    double rz = pj.r[2] - pi.r[2];

    double vx = pj.v[0] - pi.v[0];
    double vy = pj.v[1] - pi.v[1];
    double vz = pj.v[2] - pi.v[2];

    double r2   = rx*rx + ry*ry + rz*rz + ns->e2;
    //double r2   = rx*rx + ry*ry + rz*rz;
    double rinv = 1.0/sqrt(r2);

    double rv = rx*vx + ry*vy + rz*vz;

    // Add velocity
    if (sqrt(rx*rx + ry*ry + rz*rz) < hi)
    {
        nb_list[i][fi.nb] = j;
        fi.nb++;
    }

    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * pj.m;
    double mr5inv = r5inv * pj.m;

    fi.a[0] += (rx * mr3inv);
    fi.a[1] += (ry * mr3inv);
    fi.a[2] += (rz * mr3inv);

    fi.a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    fi.a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    fi.a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}



void Hermite4CPU::init_acc_jrk(Predictor *p, Forces *f, double *r_sphere)
{
    int i,j;
    #pragma omp parallel for private(j)
    for (i = 0; i < ns->n; i++)
    {
        Forces ff = {0};
        Predictor pi = p[i];
        double r_spherei = r_sphere[i];

        for (j = 0; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(pi, p[j], ff, i, j, r_spherei);
        }

        // Write on the Forces array
        f[i] = ff;
    }
}

void Hermite4CPU::update_acc_jrk(int nact, int *move, double *r_sphere,
                                 Predictor *p, Forces *f)
{
    ns->gtime.update_ini = omp_get_wtime();
    int i, j, k;
    #pragma omp parallel for private(i,j)
    for (k = 0; k < nact; k++)
    {
        i = move[k];
        // Initialization of all the struct to zero
        Forces ff = {0};
        Predictor pi = p[i];
        double r_spherei = r_sphere[i];

        #pragma omp parallel for
        for (j = 0; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(pi, p[j], ff, i, j, r_spherei);
        }

        // Write on the Forces array
        f[i] = ff;
    }
    ns->gtime.update_end += omp_get_wtime() - ns->gtime.update_ini;
}

void Hermite4CPU::predicted_pos_vel(double itime, double *t, double4 *r, double4 *v,
                                    Forces *f, Predictor *p)
{

    ns->gtime.prediction_ini = omp_get_wtime();
    for (int i = 0; i < ns->n; i++)
    {
        double dt  = itime - t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        double dt2c = dt2/2.0;
        double dt3c = dt3/6.0;

        // Temporal predictor variable
        Predictor pp;

        // Access particle information only once
        Forces ff  = f[i];
        double4 vv = v[i];
        double4 rr = r[i];

        pp.r[0] = (dt3c * ff.a1[0]) + (dt2c * ff.a[0]) + (dt * vv.x) + rr.x;
        pp.r[1] = (dt3c * ff.a1[1]) + (dt2c * ff.a[1]) + (dt * vv.y) + rr.y;
        pp.r[2] = (dt3c * ff.a1[2]) + (dt2c * ff.a[2]) + (dt * vv.z) + rr.z;

        pp.v[0] = (dt2c * ff.a1[0]) + (dt * ff.a[0]) + vv.x;
        pp.v[1] = (dt2c * ff.a1[1]) + (dt * ff.a[1]) + vv.y;
        pp.v[2] = (dt2c * ff.a1[2]) + (dt * ff.a[2]) + vv.z;

        pp.m = rr.w;

        // Update real predictor
        p[i] = pp;

    }
    ns->gtime.prediction_end += omp_get_wtime() - ns->gtime.prediction_ini;
}

void Hermite4CPU::correction_pos_vel(double itime, int nact, double *dt, double *t,
                                     int *move, Predictor *p, Forces *f, Forces *old,
                                     double3 *a2, double3 *a3,
                                     double4 *r, double4 *v)
{
    ns->gtime.correction_ini = omp_get_wtime();
    for (int k = 0; k < nact; k++)
    {
        int i = move[k];

        double dt1 = dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        double dt3c = dt3/6.0;
        double dt4c = dt4/24.0;
        double dt5c = dt5/120.0;

        double dt2inv = 1.0/dt2;
        double dt3inv = 1.0/dt3;

        // Load of the data that we will use only once, and not in every line
        Forces fi = f[i];
        Forces oldi = old[i];
        Predictor pi = p[i];
        double3 a2i;
        double3 a3i;

        // Acceleration 2nd derivate
        a2i.x = (-6 * (oldi.a[0] - fi.a[0]) - dt1 * (4 * oldi.a1[0] + 2 * fi.a1[0])) * dt2inv;
        a2i.y = (-6 * (oldi.a[1] - fi.a[1]) - dt1 * (4 * oldi.a1[1] + 2 * fi.a1[1])) * dt2inv;
        a2i.z = (-6 * (oldi.a[2] - fi.a[2]) - dt1 * (4 * oldi.a1[2] + 2 * fi.a1[2])) * dt2inv;

        // Acceleration 3rd derivate
        a3i.x = (12 * (oldi.a[0] - fi.a[0] ) + 6 * dt1 * (oldi.a1[0] + fi.a1[0])) * dt3inv;
        a3i.y = (12 * (oldi.a[1] - fi.a[1] ) + 6 * dt1 * (oldi.a1[1] + fi.a1[1])) * dt3inv;
        a3i.z = (12 * (oldi.a[2] - fi.a[2] ) + 6 * dt1 * (oldi.a1[2] + fi.a1[2])) * dt3inv;

        // Correcting position
        r[i].x = pi.r[0] + dt4c * a2i.x + dt5c * a3i.x;
        r[i].y = pi.r[1] + dt4c * a2i.y + dt5c * a3i.y;
        r[i].z = pi.r[2] + dt4c * a2i.z + dt5c * a3i.z;

        // Correcting velocity
        v[i].x = pi.v[0] + dt3c * a2i.x + dt4c * a3i.x;
        v[i].y = pi.v[1] + dt3c * a2i.y + dt4c * a3i.y;
        v[i].z = pi.v[2] + dt3c * a2i.z + dt4c * a3i.z;

        // Write on the arrays
        a2[i] = a2i;
        a3[i] = a3i;

        t[i] = itime;
        double normal_dt  = nu->get_timestep_normal(i, ns->eta);
        normal_dt = nu->normalize_dt(normal_dt, dt[i], t[i], i);

        // Ghost particles never will update its time step
        // and it will remain as the minimum possible
        //if (!ghosts[i])
            dt[i] = normal_dt;
        //else
        //    dt[i] = D_TIME_MIN;

        // TODO
        // It's required to restart the neighbours always?
        f[i].nb = 0;

    }
    ns->gtime.correction_end += omp_get_wtime() - ns->gtime.correction_ini;
}

void Hermite4CPU::update_neighbour_radius()
{
    #pragma omp parallel for
    for (int i = 0; i < ns->n; i++)
    {
        double mass = ns->h_r[i].w + ns->m_g;
        ns->h_r_sphere[i] = 5 * sqrt((ns->n * mass) * 0.5) * ns->r_cl;
    }
}

void Hermite4CPU::print_nb(double itime, int **nb_list, Forces *f, int n, Predictor *p, double *r_sphere)
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

    ns->h_dt[id] = D_TIME_MIN; // Minimum time step for new ghost particle
    // TODO: And what about the time of the particle? at the moment, we keep the same

    // Setting a zero mass to the second member.
    ns->h_r[ms.parts[1].id].w = 0.0;
    // TODO: Maybe add a blacklist to the method "find_particles_to_move"
    // using all the particles that have no mass.
}

void Hermite4CPU::multiple_systems_integration(std::vector<MultipleSystem> ms, double ITIME)
{

    cout << "Starting binary evolution for " << (int)ms.size() << " sytems\n";

    for (int i = 0; i < (int)ms.size(); i++)
    {
        // ITIME = Final time
        // CTIME = Current time
        double CTIME = 1e6;
        ms[i].next_ctime();
        ms[i].get_orbital_elements();
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
