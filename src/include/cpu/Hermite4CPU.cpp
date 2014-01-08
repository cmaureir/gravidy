#include "Hermite4CPU.hpp"

Hermite4CPU::~Hermite4CPU()
{
    // None
}

void Hermite4CPU::force_calculation(int i, int j)
{
    double rx = ns->h_p[j].r[0] - ns->h_p[i].r[0];
    double ry = ns->h_p[j].r[1] - ns->h_p[i].r[1];
    double rz = ns->h_p[j].r[2] - ns->h_p[i].r[2];

    double vx = ns->h_p[j].v[0] - ns->h_p[i].v[0];
    double vy = ns->h_p[j].v[1] - ns->h_p[i].v[1];
    double vz = ns->h_p[j].v[2] - ns->h_p[i].v[2];

    double r2     = rx*rx + ry*ry + rz*rz + ns->e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * ns->h_p[j].m;
    double mr5inv = r5inv * ns->h_p[j].m;

    double rv = rx*vx + ry*vy + rz*vz;

    ns->h_f[i].a[0] += (rx * mr3inv);
    ns->h_f[i].a[1] += (ry * mr3inv);
    ns->h_f[i].a[2] += (rz * mr3inv);

    ns->h_f[i].a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    ns->h_f[i].a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    ns->h_f[i].a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void Hermite4CPU::init_acc_jrk()
{
    int i,j;
    #pragma omp parallel for private(j) schedule(dynamic, 24)
    for (i = 0; i < ns->n; i++)
    {
        for (j = 0; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j);
        }
    }
}

void Hermite4CPU::update_acc_jrk(int nact)
{
    ns->gtime.update_ini = omp_get_wtime();
    int i, j, k;
    #pragma omp parallel for private(i,j)
    for (k = 0; k < nact; k++)
    {
        i = ns->h_move[k];
        ns->h_f[i].a[0]  = 0.0;
        ns->h_f[i].a[1]  = 0.0;
        ns->h_f[i].a[2]  = 0.0;
        ns->h_f[i].a1[0] = 0.0;
        ns->h_f[i].a1[1] = 0.0;
        ns->h_f[i].a1[2] = 0.0;

        #pragma omp parallel for
        for (j = 0; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j);
        }
    }
    ns->gtime.update_end += omp_get_wtime() - ns->gtime.update_ini;
}

void Hermite4CPU::predicted_pos_vel(double ITIME)
{

    ns->gtime.prediction_ini = omp_get_wtime();
    for (int i = 0; i < ns->n; i++)
    {
        double dt  = ITIME - ns->h_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        ns->h_p[i].r[0] = (dt3/6 * ns->h_f[i].a1[0]) + (dt2/2 * ns->h_f[i].a[0]) + (dt * ns->h_v[i].x) + ns->h_r[i].x;
        ns->h_p[i].r[1] = (dt3/6 * ns->h_f[i].a1[1]) + (dt2/2 * ns->h_f[i].a[1]) + (dt * ns->h_v[i].y) + ns->h_r[i].y;
        ns->h_p[i].r[2] = (dt3/6 * ns->h_f[i].a1[2]) + (dt2/2 * ns->h_f[i].a[2]) + (dt * ns->h_v[i].z) + ns->h_r[i].z;

        ns->h_p[i].v[0] = (dt2/2 * ns->h_f[i].a1[0]) + (dt * ns->h_f[i].a[0]) + ns->h_v[i].x;
        ns->h_p[i].v[1] = (dt2/2 * ns->h_f[i].a1[1]) + (dt * ns->h_f[i].a[1]) + ns->h_v[i].y;
        ns->h_p[i].v[2] = (dt2/2 * ns->h_f[i].a1[2]) + (dt * ns->h_f[i].a[2]) + ns->h_v[i].z;

        ns->h_p[i].m = ns->h_r[i].w;

    }
    ns->gtime.prediction_end += omp_get_wtime() - ns->gtime.prediction_ini;
}

void Hermite4CPU::correction_pos_vel(double ITIME, int nact)
{
    ns->gtime.correction_ini = omp_get_wtime();
    for (int k = 0; k < nact; k++)
    {
        int i = ns->h_move[k];

        double dt1 = ns->h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        ns->h_a2[i].x = (-6 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) - dt1 * (4 * ns->h_old[i].a1[0] + 2 * ns->h_f[i].a1[0]) ) / dt2;
        ns->h_a2[i].y = (-6 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) - dt1 * (4 * ns->h_old[i].a1[1] + 2 * ns->h_f[i].a1[1]) ) / dt2;
        ns->h_a2[i].z = (-6 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) - dt1 * (4 * ns->h_old[i].a1[2] + 2 * ns->h_f[i].a1[2]) ) / dt2;

        // Acceleration 3rd derivate
        ns->h_a3[i].x = (12 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) + 6 * dt1 * (ns->h_old[i].a1[0] + ns->h_f[i].a1[0]) ) / dt3;
        ns->h_a3[i].y = (12 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) + 6 * dt1 * (ns->h_old[i].a1[1] + ns->h_f[i].a1[1]) ) / dt3;
        ns->h_a3[i].z = (12 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) + 6 * dt1 * (ns->h_old[i].a1[2] + ns->h_f[i].a1[2]) ) / dt3;

        // Correcting position
        ns->h_r[i].x = ns->h_p[i].r[0] + (dt4/24)*ns->h_a2[i].x + (dt5/120)*ns->h_a3[i].x;
        ns->h_r[i].y = ns->h_p[i].r[1] + (dt4/24)*ns->h_a2[i].y + (dt5/120)*ns->h_a3[i].y;
        ns->h_r[i].z = ns->h_p[i].r[2] + (dt4/24)*ns->h_a2[i].z + (dt5/120)*ns->h_a3[i].z;

        // Correcting velocity
        ns->h_v[i].x = ns->h_p[i].v[0] + (dt3/6)*ns->h_a2[i].x + (dt4/24)*ns->h_a3[i].x;
        ns->h_v[i].y = ns->h_p[i].v[1] + (dt3/6)*ns->h_a2[i].y + (dt4/24)*ns->h_a3[i].y;
        ns->h_v[i].z = ns->h_p[i].v[2] + (dt3/6)*ns->h_a2[i].z + (dt4/24)*ns->h_a3[i].z;


        ns->h_t[i] = ITIME;
        double normal_dt  = nu->get_timestep_normal(i);
        normal_dt = nu->normalize_dt(normal_dt, ns->h_dt[i], ns->h_t[i], i);
        ns->h_dt[i] = normal_dt;

    }
    ns->gtime.correction_end += omp_get_wtime() - ns->gtime.correction_ini;
}



void Hermite4CPU::integration()
{

    ns->gtime.integration_ini = omp_get_wtime();

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;


    int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads - 1);

    init_acc_jrk();
    init_dt(ATIME);

    ns->en.ini = nu->get_energy();   // Initial calculation of the energy of the system
    ns->en.tmp = ns->en.ini;

    ns->hmr_time = nu->get_half_mass_relaxation_time();
    ns->cr_time  = nu->get_crossing_time();

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
