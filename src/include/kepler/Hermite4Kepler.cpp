#include "Hermite4Kepler.hpp"

/* @desc Constructor which uses the parent constructor (Hermite4).
 *       The additional steps in this method are the allocation
 *       and initialization of the data that it is only needed
 *       for the Keplerian treatment.
 */
Hermite4Kepler::Hermite4Kepler(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
: Hermite4(ns, logger, nu)
{
    alloc_arrays_host_kepler();
    init_data_bh();
}

/* @desc Destructor in charge of freeing the memory of the data
 *       used in the Keplerian treatment
 */
Hermite4Kepler::~Hermite4Kepler()
{
    free_arrays_host_kepler();
}

void Hermite4Kepler::alloc_arrays_host_kepler()
{

    size_t ff_size = ns->n * sizeof(Forces);
    size_t d4_size = ns->n * sizeof(double4);
    size_t d1_size = ns->n * sizeof(double);

    ns->h_fbh    = new Forces[ff_size];
    ns->h_a2bh   = new double4[d4_size];
    ns->h_a3bh   = new double4[d4_size];
    ns->h_dt_old = new double[d1_size];
}

void Hermite4Kepler::free_arrays_host_kepler()
{
    delete ns->h_fbh;
    delete ns->h_a2bh;
    delete ns->h_a2bh;
    delete ns->h_dt_old;
}

void Hermite4Kepler::init_data_bh()
{
    double4 empty = {0.0, 0.0, 0.0, 0.0};
    int i;
    #pragma omp parallel for schedule(dynamic, 24)
    for (i = 1; i < ns->n; i++) {
        ns->h_fbh[i].a[0]  = 0.0;
        ns->h_fbh[i].a[1]  = 0.0;
        ns->h_fbh[i].a[2]  = 0.0;

        ns->h_fbh[i].a1[0] = 0.0;
        ns->h_fbh[i].a1[1] = 0.0;
        ns->h_fbh[i].a1[2] = 0.0;

        ns->h_a2bh[i] = empty;
        ns->h_a3bh[i] = empty;
        ns->h_dt_old[i] = 0.0;
    }
}

void Hermite4Kepler::force_calculation(int i, int j)
{
    double rx = ns->h_p[j].r[0] - ns->h_p[i].r[0];
    double ry = ns->h_p[j].r[1] - ns->h_p[i].r[1];
    double rz = ns->h_p[j].r[2] - ns->h_p[i].r[2];

    double vx = ns->h_p[j].v[0] - ns->h_p[i].v[0];
    double vy = ns->h_p[j].v[1] - ns->h_p[i].v[1];
    double vz = ns->h_p[j].v[2] - ns->h_p[i].v[2];

    double r2     = (rx * rx + ry * ry) + (rz * rz + ns->e2);
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

/*
 * @name force_calculation_bh
 *
 * @desc Method to calculate the gravitational interaction between
 * a certain i-particle and the BH
 */
void Hermite4Kepler::force_calculation_bh(int i)
{
    double rx = ns->h_p[0].r[0] - ns->h_p[i].r[0];
    double ry = ns->h_p[0].r[1] - ns->h_p[i].r[1];
    double rz = ns->h_p[0].r[2] - ns->h_p[i].r[2];

    double vx = ns->h_p[0].v[0] - ns->h_p[i].v[0];
    double vy = ns->h_p[0].v[1] - ns->h_p[i].v[1];
    double vz = ns->h_p[0].v[2] - ns->h_p[i].v[2];

    double r2     = (rx * rx + ry * ry) + (rz * rz + ns->e2);
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * ns->h_p[0].m;
    double mr5inv = r5inv * ns->h_p[0].m;

    double rv = rx*vx + ry*vy + rz*vz;

    ns->h_fbh[i].a[0] += (rx * mr3inv);
    ns->h_fbh[i].a[1] += (ry * mr3inv);
    ns->h_fbh[i].a[2] += (rz * mr3inv);

    ns->h_fbh[i].a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    ns->h_fbh[i].a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    ns->h_fbh[i].a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void Hermite4Kepler::init_acc_jrk()
{
    int i,j;
    #pragma omp parallel for private(j) schedule(dynamic, 24)
    for (i = 1; i < ns->n; i++)
    {
        for (j = 1; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j);
        }
    }
}

void Hermite4Kepler::init_acc_jrk_bh()
{
    int i;
    #pragma omp parallel for schedule(dynamic, 24)
    for (i = 1; i < ns->n; i++)
    {
        force_calculation_bh(i);
    }
}

void Hermite4Kepler::update_acc_jrk(int nact)
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
        for (j = FIRST_PARTICLE; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j);
        }
    }
    ns->gtime.update_end += omp_get_wtime() - ns->gtime.update_ini;
}

void Hermite4Kepler::predicted_pos_vel(double ITIME)
{

    ns->gtime.prediction_ini = omp_get_wtime();
    for (int i = 1; i < ns->n; i++)
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

void Hermite4Kepler::correction_pos_vel(double ITIME, int nact)
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


        // Save old timestep
        ns->h_dt_old[i] = ns->h_dt[i];
        // Update self timej
        ns->h_t[i] = ITIME;

        // Get new timestep
        double normal_dt  = nu->get_timestep_normal(i);
        normal_dt = nu->normalize_dt(normal_dt, ns->h_dt[i], ns->h_t[i], i);
        ns->h_dt[i] = normal_dt;

    }
    ns->gtime.correction_end += omp_get_wtime() - ns->gtime.correction_ini;
}

void Hermite4Kepler::integration()
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
    init_acc_jrk_bh();
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

        move_keplerian_orbit(ITIME, nact);
        predicted_pos_vel(ITIME);

        update_acc_jrk(nact);

        correction_pos_vel(ITIME, nact);

        // Update the amount of interactions counter
        interactions += nact * ns->n;

        // Find the next integration time
        next_integration_time(ATIME);

        /****************************/
        //std::cout << nact << " | ";
        //for (int ii = 0; ii < nact; ii++) {
        //    std::cout << ns->h_move[ii] << " ";
        //}
        //std::cout << std::endl;
        /****************************/

        if(std::ceil(ITIME) == ITIME)
        {
            assert(nact == ns->n - FIRST_PARTICLE);
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

// -------------------------------------------------------------------------------

/*
 * @fn move_keplerian_orbit()
 *
 * @brief
 *  Perform the calculation of the predicted position
 *    and velocity using the Kepler's equation.
 *
 */
void Hermite4Kepler::move_keplerian_orbit(double ITIME, int nact)
{
    for (int i = 0; i < nact; i++)
    {
        int k = ns->h_move[0];
        double dt = ITIME - ns->h_t[k];
        double time = 0.0;

        pkepler.r[0] = ns->h_r[k].x - ns->h_r[0].x;
        pkepler.r[1] = ns->h_r[k].y - ns->h_r[0].y;
        pkepler.r[2] = ns->h_r[k].z - ns->h_r[0].z;

        pkepler.v[0] = ns->h_v[k].x - ns->h_v[0].x;
        pkepler.v[1] = ns->h_v[k].y - ns->h_v[0].y;
        pkepler.v[2] = ns->h_v[k].z - ns->h_v[0].z;
        pkepler.m    = ns->h_r[0].w;

        for (time = ns->h_t[k]; time < ITIME; time += dt)
        {
            kepler_move(k, time);
        }

        ns->h_p[k] = pkepler;
    }
}
void Hermite4Kepler::calculate_orbital_elements()
{
    double b_mag;
    double mu;


    mu = G * pkepler.m;

    r_mag = sqrt((pkepler.r[0]) * (pkepler.r[0]) +
                 (pkepler.r[1]) * (pkepler.r[1]) +
                 (pkepler.r[2]) * (pkepler.r[2]));

    // Angular momentum
    // j = r x v
    oe.j.x = (pkepler.r[1]) * (pkepler.v[2]) - (pkepler.r[2]) * (pkepler.v[1]);
    oe.j.y = (pkepler.r[2]) * (pkepler.v[0]) - (pkepler.r[0]) * (pkepler.v[2]);
    oe.j.z = (pkepler.r[0]) * (pkepler.v[1]) - (pkepler.r[1]) * (pkepler.v[0]);

    // Runge-Lenz-vector
    // e = { (v x j) / (G * m) }  - { r / |r| }
    oe.e.x = (pkepler.v[1]) * (oe.j.z) - (pkepler.v[2]) * (oe.j.y);
    oe.e.y = (pkepler.v[2]) * (oe.j.x) - (pkepler.v[0]) * (oe.j.z);
    oe.e.z = (pkepler.v[0]) * (oe.j.y) - (pkepler.v[1]) * (oe.j.x);

    // G, omitted due its value of 1
    oe.e.x = oe.e.x/mu - pkepler.r[0]/r_mag;
    oe.e.y = oe.e.y/mu - pkepler.r[1]/r_mag;
    oe.e.z = oe.e.z/mu - pkepler.r[2]/r_mag;

    // Eccentricity
    oe.ecc = sqrt(oe.e.x * oe.e.x +
                  oe.e.y * oe.e.y +
                  oe.e.z * oe.e.z);

    // Semi-major axis
    // a = ( j * j ) / (G * m * | 1 - ecc^2 | )
    oe.a = (oe.j.x * oe.j.x +
            oe.j.y * oe.j.y +
            oe.j.z * oe.j.z);
    oe.a /= (mu * std::abs(1 - oe.ecc * oe.ecc));

    // Frequency of the orbit
    // w = sqrt(( G * m )/ a^3 )
    oe.w = sqrt( mu / (oe.a * oe.a * oe.a));

    // Semi-major vector
    oe.a_vec.x = oe.a * oe.e.x/oe.ecc;
    oe.a_vec.y = oe.a * oe.e.y/oe.ecc;
    oe.a_vec.z = oe.a * oe.e.z/oe.ecc;

    // Semi-minor axis
    oe.b = oe.a * sqrt(std::abs(1 - oe.ecc * oe.ecc));

    // Semi-minor vector
    oe.b_vec.x = (oe.j.y * oe.e.z - oe.j.z * oe.e.y);
    oe.b_vec.y = (oe.j.z * oe.e.x - oe.j.x * oe.e.z);
    oe.b_vec.z = (oe.j.x * oe.e.y - oe.j.y * oe.e.x);

    b_mag = sqrt(oe.b_vec.x * oe.b_vec.x +
                 oe.b_vec.y * oe.b_vec.y +
                 oe.b_vec.z * oe.b_vec.z);

    oe.b_vec.x *= oe.b/b_mag;
    oe.b_vec.y *= oe.b/b_mag;
    oe.b_vec.z *= oe.b/b_mag;

    oe.e_anomaly = 0.0; // The value will depend of the orbit type
    oe.m_anomaly = 0.0; // The value will depend of the orbit type
}

void Hermite4Kepler::print_orbital_elements()
{
    int ll = 30;
    std::cout.precision(15);
    std::cout << std::setw(ll) << "Orbital Elements:" << std::endl;

    std::cout << std::setw(ll) << "Angular Momentum vector (";
    std::cout << oe.j.x << ", " << oe.j.y << ", " << oe.j.z;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Runge-Lenz vector (";
    std::cout << oe.e.x << ", " << oe.e.y << ", " << oe.e.z;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Semi-major axis (";
    std::cout << oe.a;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Semi-minor axis (";
    std::cout << oe.b;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Frequency of the orbit (";
    std::cout << oe.w;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Semi-major axis vector (";
    std::cout << oe.a_vec.x << ", " << oe.a_vec.y << ", " << oe.a_vec.z;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Semi-minor axis vector (";
    std::cout << oe.b_vec.x << ", " << oe.b_vec.y << ", " << oe.b_vec.z;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Eccentricity (";
    std::cout << oe.ecc;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Eccentric anomaly (";
    std::cout << oe.e_anomaly;
    std::cout << ")" << std::endl;

    std::cout << std::setw(ll) << "Mean anomaly (";
    std::cout << oe.m_anomaly;
    std::cout << ")" << std::endl;
}

Predictor Hermite4Kepler::get_elliptical_pos_vel(int i, double dt)
{
    double r_const;
    double v_const;
    double cos_e;
    double sin_e;
    double ecc_const;
    double cos_const;

    //if (dt == 0)
    if (dt == ns->h_t[i])
    {
        // Calculate Cosine of Eccentric Anomaly
        oe.e_anomaly0 = (oe.a - r_mag) / (oe.ecc * oe.a);
        //std::cout << "e_anomaly(0): " << oe.e_anomaly << std::endl;

        // Fixing Cosine argument
        if(oe.e_anomaly0 >= 1.0)
            oe.e_anomaly0 = 0.0;
        else if(oe.e_anomaly0 <= -1.0)
            oe.e_anomaly0 = M_PI;
        else
            oe.e_anomaly0 = acos(oe.e_anomaly0);

        oe.m_anomaly0 = (oe.e_anomaly0 - oe.ecc*sin(oe.e_anomaly0));
    }

    // Calculate Mean Anomaly
    oe.m_anomaly = oe.m_anomaly0 + dt * oe.w;

    // Adjusting M anomaly to be < 2 * pi
    if(oe.m_anomaly >= 2 * M_PI)
        oe.m_anomaly = fmod(oe.m_anomaly, 2 * M_PI);

    // Solving the Kepler equation for elliptical orbits
    oe.e_anomaly = solve_kepler(oe.m_anomaly, oe.ecc);

    cos_e = cos(oe.e_anomaly);
    sin_e = sin(oe.e_anomaly);

    r_const = cos_e - oe.ecc;
    v_const = oe.w / (1.0 - oe.ecc * cos_e);

    // Based on Loeckmann and Baumgardt simulations criteria
    // This work better with 0.99 < e < 1 and |E| < 1e-3
    if(oe.ecc > 0.99)
    {
        double e_tmp = (oe.e_anomaly > 2.0 * M_PI - 1e-3) ? oe.e_anomaly - 2.0 * M_PI : oe.e_anomaly;
        if(e_tmp < 1e-3)
        {
            e_tmp *= e_tmp;
            double j_mag = oe.j.x * oe.j.x +
                           oe.j.y * oe.j.y +
                           oe.j.z * oe.j.z;
            ecc_const = j_mag/(ns->h_r[0].w * oe.a * (1 + oe.ecc));
            cos_const = -0.5 * e_tmp * (1 - e_tmp / 12.0 * (1 - e_tmp / 30.0));

            r_const = ecc_const + cos_const;
            v_const = oe.w / (ecc_const - oe.ecc * cos_const);

        }
    }

    Predictor ppv;
    // New position
    ppv.r[0] =   oe.a_vec.x * r_const + oe.b_vec.x * sin_e;
    ppv.r[1] =   oe.a_vec.y * r_const + oe.b_vec.y * sin_e;
    ppv.r[2] =   oe.a_vec.z * r_const + oe.b_vec.z * sin_e;

    // New velocity
    ppv.v[0] = (-oe.a_vec.x * sin_e + oe.b_vec.x * cos_e) * v_const;
    ppv.v[1] = (-oe.a_vec.y * sin_e + oe.b_vec.y * cos_e) * v_const;
    ppv.v[2] = (-oe.a_vec.z * sin_e + oe.b_vec.z * cos_e) * v_const;

    return ppv;
}

Predictor Hermite4Kepler::get_hyperbolic_pos_vel(int i, double dt)
{

    //if (dt == 0)
    if (dt == ns->h_t[i])
    {
        rdotv = (pkepler.r[0] * pkepler.v[0] +
                 pkepler.r[1] * pkepler.v[1] +
                 pkepler.r[2] * pkepler.v[2]);
        // calculate eccentric anomaly e at t+dt
        oe.e_anomaly0 = (oe.a + r_mag) / (oe.ecc * oe.a);

        if(oe.e_anomaly0 < 1.0)
            oe.e_anomaly0 = 0.0;
        else if(rdotv < 0)
            oe.e_anomaly0 = -acosh(oe.e_anomaly0);
        else
            oe.e_anomaly0 = acosh(oe.e_anomaly0);

        oe.m_anomaly0 = oe.ecc * sinh(oe.e_anomaly0) - oe.e_anomaly0;
    }

    oe.m_anomaly = oe.m_anomaly0 + dt * oe.w;

    // Solving Kepler's equation for Hyperbolic/Parabolic orbits
    oe.e_anomaly = kepler(oe.ecc, oe.m_anomaly);

    double cos_e = cosh(oe.e_anomaly);
    double sin_e = sinh(oe.e_anomaly);
    double v_const = oe.w / (oe.ecc * cos_e - 1.);

    Predictor ppv;
    // New position
    ppv.r[0] = oe.a_vec.x * (oe.ecc - cos_e)  + oe.b_vec.x * sin_e;
    ppv.r[1] = oe.a_vec.y * (oe.ecc - cos_e)  + oe.b_vec.y * sin_e;
    ppv.r[2] = oe.a_vec.z * (oe.ecc - cos_e)  + oe.b_vec.z * sin_e;

    // New velocity
    ppv.v[0] = (-oe.a_vec.x * sin_e + oe.b_vec.x * cos_e) * v_const;
    ppv.v[1] = (-oe.a_vec.y * sin_e + oe.b_vec.y * cos_e) * v_const;
    ppv.v[2] = (-oe.a_vec.z * sin_e + oe.b_vec.z * cos_e) * v_const;

    return ppv;
}

void Hermite4Kepler::kepler_move(int i, double dt)
{

    // Calculating some orbital elements of the interaction between an i-particle
    // and a BH */
    calculate_orbital_elements();

    //print_orbital_elements();
    //getchar();
    Predictor new_rv;

    if(oe.ecc < 1)
    {
        new_rv = get_elliptical_pos_vel(i, dt);
    }
    else
    {
        new_rv = get_hyperbolic_pos_vel(i, dt);
    }

    // Hacer algo con pkepler y el valor que obtuvimos de new_pos_vel
    //
    // End of the orbit-type treatment

    //double new_r_mag = sqrt((new_rv.r[0]) * (new_rv.r[0]) +
    //                        (new_rv.r[1]) * (new_rv.r[1]) +
    //                        (new_rv.r[2]) * (new_rv.r[2]));
    //double new_r2 = new_r_mag * new_r_mag;

    //double new_v_mag = sqrt((new_rv.v[0]) * (new_rv.v[0]) +
    //                        (new_rv.v[1]) * (new_rv.v[1]) +
    //                        (new_rv.v[2]) * (new_rv.v[2]));
    //double new_v2 = new_v_mag * new_v_mag;

    //double vr = ((new_rv.v[0]) * (new_rv.r[0]) +
    //             (new_rv.v[1]) * (new_rv.r[1]) +
    //             (new_rv.v[2]) * (new_rv.r[2]));

    //double j_mag = sqrt(oe.j.x * oe.j.x +
    //                    oe.j.y * oe.j.y +
    //                    oe.j.z * oe.j.z);

    //double rv2 = (new_rv.r[0] * new_rv.v[0] +
    //              new_rv.r[1] * new_rv.v[1] +
    //              new_rv.r[2] * new_rv.v[2]);

    // get |v| from j = r x v
    // Esta mierda no se usa...
    //new_v_mag = j_mag / (new_r_mag * new_v_mag * sin(acos(rv2/(new_r_mag * new_v_mag))));

    // total motion relative to fix central mass
    //pkepler.r[0] = ns->h_r[0].x + (new_rv.r[0]);
    //pkepler.r[1] = ns->h_r[0].y + (new_rv.r[1]);
    //pkepler.r[2] = ns->h_r[0].z + (new_rv.r[2]);

    //pkepler.v[0] = ns->h_v[0].x + (new_rv.v[0]);
    //pkepler.v[0] = ns->h_v[0].y + (new_rv.v[1]);
    //pkepler.v[0] = ns->h_v[0].z + (new_rv.v[2]);
    pkepler.r[0] =  new_rv.r[0];
    pkepler.r[1] =  new_rv.r[1];
    pkepler.r[2] =  new_rv.r[2];

    pkepler.v[0] =  new_rv.v[0];
    pkepler.v[1] =  new_rv.v[1];
    pkepler.v[2] =  new_rv.v[2];

    ///*
    // * Force contribution of the central mass on a particle
    // * (Acceleration and it 1st, 2nd and 3rd derivatives
    // */

    //double r2inv = 1/(r2);
    //double mr3inv = - ns->h_r[0].w * r2inv * sqrt(r2inv);

    //// Acceleration (a)
    //// a = -G * ns->h_r[0].w * \vec{r} / r^{3}

    //double a0x =  (new_rv.r[0]) * mr3inv;
    //double a0y =  (new_rv.r[1]) * mr3inv;
    //double a0z =  (new_rv.r[2]) * mr3inv;


    //double a1x = mr3inv * ( (new_rv.v[0]) - 3 * r2inv * vr * (new_rv.r[0]));
    //double a1y = mr3inv * ( (new_rv.v[1]) - 3 * r2inv * vr * (new_rv.r[1]));
    //double a1z = mr3inv * ( (new_rv.v[2]) - 3 * r2inv * vr * (new_rv.r[2]));

    //double ra0 = ((new_rv.r[0])*a0x + (new_rv.r[1])*a0y + (new_rv.r[2])*a0z);

    //double a2x = mr3inv * (a0x - 3.0 * r2inv * (vr * ( 2.0 * (new_rv.v[0]) - 5.0 * vr * (new_rv.r[0]) * r2inv)) + (v2 + ra0) * (new_rv.r[0]));
    //double a2y = mr3inv * (a0y - 3.0 * r2inv * (vr * ( 2.0 * (new_rv.v[1]) - 5.0 * vr * (new_rv.r[1]) * r2inv)) + (v2 + ra0) * (new_rv.r[1]));
    //double a2z = mr3inv * (a0z - 3.0 * r2inv * (vr * ( 2.0 * (new_rv.v[2]) - 5.0 * vr * (new_rv.r[2]) * r2inv)) + (v2 + ra0) * (new_rv.r[2]));

    //double va = ((new_rv.v[0])*a0x + (new_rv.v[1])*a0y + (new_rv.v[2])*a0z);
    //double ra1 = ((new_rv.r[0])*a1x + (new_rv.r[1])*a1y + (new_rv.r[2])*a1z);

    //double a3x = mr3inv * ( a1x - 3.0 * r2inv * ( 3.0 * vr * a0x + 3.0
    //            * (v2 + ra0) * ((new_rv.v[0]) - 5.0 * vr * r2inv * (new_rv.r[0]))
    //            + (3.0 * va + ra1) * (new_rv.r[0]) + (vr * vr * r2inv)
    //            * (-15.0 * (new_rv.v[0]) + 35.0 * vr * r2inv * (new_rv.r[0]))));

    //double a3y = mr3inv * ( a1y - 3.0 * r2inv * ( 3.0 * vr * a0y + 3.0
    //            * (v2 + ra0) * ((new_rv.v[1]) - 5.0 * vr * r2inv * (new_rv.r[1]))
    //            + (3.0 * va + ra1) * (new_rv.r[1]) + (vr * vr * r2inv)
    //            * (-15.0 * (new_rv.v[1]) + 35.0 * vr * r2inv * (new_rv.r[1]))));

    //double a3z = mr3inv * ( a1z - 3.0 * r2inv * ( 3.0 * vr * a0z + 3.0
    //            * (v2 + ra0) * ((new_rv.v[2]) - 5.0 * vr * r2inv * (new_rv.r[2]))
    //            + (3.0 * va + ra1) * (new_rv.r[2]) + (vr * vr * r2inv)
    //            * (-15.0 * (new_rv.v[2]) + 35.0 * vr * r2inv * (new_rv.r[2]))));

    //h_fbh[i].a[0] = a0x;
    //h_fbh[i].a[1] = a0y;
    //h_fbh[i].a[2] = a0z;

    //h_fbh[i].a1[0] = a1x;
    //h_fbh[i].a1[1] = a1y;
    //h_fbh[i].a1[2] = a1z;

    //h_a2bh[i].x = a2x;
    //h_a2bh[i].y = a2y;
    //h_a2bh[i].z = a2z;

    //h_a3bh[i].x = a3x;
    //h_a3bh[i].y = a3y;
    //h_a3bh[i].z = a3z;

}

double Hermite4Kepler::solve_kepler(double m_anomaly, double ecc)
{
    // Solving Kepler's equation for an elliptical orbit
    double e_new = ecc > 0.8 ? M_PI : m_anomaly;
    double d = 0;
    int e_iter = 0;

    while(fabs(d) > DEL_E)
    {
        d = e_new - ecc * sin(e_new) - m_anomaly;
        if(e_iter-1 >= KEPLER_ITE)
            break;

        e_new -= d / (1.0 - ecc * cos(e_new));
        e_iter++;
    }
    return e_new;
}

/*
 * Following function taken from
 * http://www.projectpluto.com/kepler.htm
 * (Adapted to solve only hyperbolic orbits.)
 */
double Hermite4Kepler::kepler(const double ecc, double mean_anom)
{
    double curr, err;
    int is_negative = 0, n_iter = 0;

    if(!mean_anom)
        return(0.);

    if(ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
    {
        curr = atan2(sin(mean_anom), cos(mean_anom) - ecc);
        /* one correction step,  and we're done */
        err = curr - ecc * sin(curr) - mean_anom;
        curr -= err / (1. - ecc * cos(curr));
        //printf("(1) %.10f\n", curr);
        return(curr);
    }

    if(mean_anom < 0.)
    {
        mean_anom = -mean_anom;
        is_negative = 1;
    }

    curr = mean_anom;
    if((ecc > .8 && mean_anom < M_PI / 3.) || ecc > 1.)    /* up to 60 degrees */
    {
        double trial = mean_anom / fabs(1. - ecc);

        if(trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
        {
            if(mean_anom < M_PI)
                trial = pow(6. * mean_anom, 1.0/3.0);
            else        /* hyperbolic w/ 5th & higher-order terms predominant */
                trial = asinh(mean_anom / ecc);
        }
        curr = trial;
    }

    double thresh = DEL_E_HYP * mean_anom;
    if(ecc < 1.)
    {
        err = (curr - ecc * sin(curr)) - mean_anom;
        while(fabs(err) > thresh)
        {
            n_iter++;
            curr -= err / (1. - ecc * cos(curr));
            err = (curr - ecc * sin(curr)) - mean_anom;

            if(n_iter > KEPLER_ITE) // amended
            {
                #ifdef DEBUG_KEPLER
                printf(
                    "#### ! Aborting kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%e) ####\n",
                    KEPLER_ITE, err,
                    ecc, mean_anom, curr);
                #endif
                break;
            }
        }
    }
    else
    {
        curr = log(2.*mean_anom/ecc+1.8); // taken from Burkardt & Danby, CeMec 31 (1983), 317-328
        double curr_abs = fabs(curr);
        err = (ecc * sinh(curr) - curr) - mean_anom;
        while(fabs(err) > thresh)
        {
            n_iter++;
            if(curr_abs < .72 && ecc < 1.1)
            {
                // [e * sinh(E) - E] / E << 1, and/or e * cosh(E) - 1 << 1
                // so don't calculate it directly
                double curr2 = curr * curr;

                // relative error when omitting nth order term needs to be smaller than resolution 1.e-15:
                // .5 * E^2 > 1e15 * E^n/n!, i.e. E < (n!/2e15)^(1/(n-2))
                // n = 16: E < .72, n = 10: E < .08

                if(curr_abs > .08)
                    curr -= err / ((ecc - 1) * cosh(curr) +
                            (((((((_1_16_15 * curr2 + 1.)
                            * _1_14_13 * curr2 + 1.)
                            * _1_12_11 * curr2 + 1.)
                            * _1_10_9 * curr2 + 1.)
                            * _1_8_7 * curr2 + 1.)
                            * _1_6_5 * curr2 + 1.)
                            * _1_4_3 * curr2 + 1.)
                            * .5 *curr2);
                else
                    curr -= err / ((ecc - 1) * cosh(curr) +
                            ((((_1_10_9 * curr2 + 1.)
                            * _1_8_7 * curr2 + 1.)
                            * _1_6_5 * curr2 + 1.)
                            * _1_4_3 * curr2 + 1.)
                            * .5 *curr2);
                curr2 = curr * curr;
                curr_abs = fabs(curr);

                if(curr_abs > .08)
                    err = ((ecc - 1) * sinh(curr) +
                            (((((((_1_17_16 * curr2 + 1.)
                            * _1_15_14 * curr2 + 1.)
                            * _1_13_12 * curr2 + 1.)
                            * _1_11_10 * curr2 + 1.)
                            * _1_9_8 * curr2 + 1.)
                            * _1_7_6 * curr2 + 1.)
                            * .05 * curr2 + 1.)
                            * _1_3_2 * curr2 * curr) - mean_anom;
                else
                    err = ((ecc - 1) * sinh(curr) +
                            ((((_1_11_10 * curr2 + 1.)
                            * _1_9_8 * curr2 + 1.)
                            * _1_7_6 * curr2 + 1.)
                            * .05 * curr2 + 1.)
                            * _1_3_2 * curr2 * curr) - mean_anom;
            }
            else
            {
                curr -= err / (ecc * cosh(curr) - 1.);
                err = (ecc * sinh(curr) - curr) - mean_anom;
            }
            if(n_iter > KEPLER_ITE) // amended
            {
                #ifdef DEBUG_KEPLER
                printf(
                    "### Aborting hyperbolic kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%1.12e, sinh(E_)=%1.12e)\n",
                    KEPLER_ITE, err,
                    ecc, mean_anom, curr, sinh(curr));
                #endif

                break;
            }
        }
     }

    return(is_negative ? -curr : curr);
}
