#include "dynamics_cpu.hpp"

/*
 * @fn next_itime()
 *
 * @brief
 *  Determine the next iteration time
 *  considering the time of all the particles
 */
void next_itime(double *ATIME)
{
    // Big number to find the minimum
    *ATIME = 1.0e10;
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double time = h_t[i] + h_dt[i];
        if(time < *ATIME)
        {
            *ATIME = time;
        }
    }
}

/*
 * @fn find_particles_to_move()
 *
 * @brief
 *  Find particles wich have the same time
 *  that the integration time defined in the
 *  function 'next_time'
 */
int find_particles_to_move(double ITIME)
{
    int j = 0;
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        h_move[i] = -1;
        if (h_t[i] + h_dt[i] == ITIME)
        {
            h_move[j] = i;
            j++;
        }
    }
    return j;
}

/*
 * @fn init_dt()
 *
 * @brief
 *  Initialization of each particle time step,
 *  following the Aarseth approach.
 *
 *  After set each time step, it is necessah_ry
 *  to adjust it to the block time steps.
 */
void init_dt(double *ATIME)
{
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    double tmp_dt;
    for (int i = 0; i < n; i++)
    {
        double a2 = get_magnitude(h_a[i].x, h_a[i].y, h_a[i].z);
        double j2 = get_magnitude(h_j[i].x, h_j[i].y, h_j[i].z);
        tmp_dt = ETA_S * (a2/j2);

        // Adjusting to block timesteps
        // to the nearest-lower power of two
        int exp = (int)(std::ceil(log(tmp_dt)/log(2.0))-1);
        tmp_dt = pow(2,exp);

        if (tmp_dt < D_TIME_MIN)
            tmp_dt = D_TIME_MIN;
        else if (tmp_dt > D_TIME_MAX)
            tmp_dt = D_TIME_MAX;

        h_dt[i] = tmp_dt;
        h_t[i] = 0.0;

        // Obtaining the first integration time
        if(tmp_dt < *ATIME)
            *ATIME = tmp_dt;
    }
}

void force_calculation(int i, int j)
{
    double rx = h_p_r[j].x - h_p_r[i].x;
    double ry = h_p_r[j].y - h_p_r[i].y;
    double rz = h_p_r[j].z - h_p_r[i].z;

    double vx = h_p_v[j].x - h_p_v[i].x;
    double vy = h_p_v[j].y - h_p_v[i].y;
    double vz = h_p_v[j].z - h_p_v[i].z;

    double r2 = rx*rx + ry*ry + rz*rz + E2;
    double rinv = 1/sqrt(r2);
    double r2inv = rinv  * rinv;
    double r3inv = r2inv * rinv;
    double r5inv = r2inv * r3inv;
    double mr3inv = r3inv * h_m[j];
    double mr5inv = r5inv * h_m[j];

    double rv = rx*vx + ry*vy + rz*vz;

    h_a[i].x += (rx * mr3inv);
    h_a[i].y += (ry * mr3inv);
    h_a[i].z += (rz * mr3inv);

    h_j[i].x += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    h_j[i].y += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    h_j[i].z += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void init_acc_jrk()
{
    int j;
//    #pragma omp parallel for private(j)
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i,j);
        }

    }
}

/*
 * @fn update_acc_jrk()
 *
 * @brief
 *  Iteration calculation of the acceleration and the jrk.
 *
 * @todo
 *  Only perform the calculation of the particles with the t + dt = ITIME
 *
 */
void update_acc_jrk(int total)
{
    int i, j;
    //#pragma omp parallel for private(i,j)
    for (int k = 0; k < total; k++)
    {
        i = h_move[k];
        // Cleaning acceleration and jrk
        h_a[i].x = 0.0;
        h_a[i].y = 0.0;
        h_a[i].z = 0.0;
        h_j[i].x = 0.0;
        h_j[i].y = 0.0;
        h_j[i].z = 0.0;

        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i,j);
        }

    //    printf("%d %.10f %.10f %.10f - %.10f %.10f %.10f\n",
    //            i, h_a[i].x, h_a[i].y, h_a[i].z,
    //            h_j[i].x, h_j[i].y, h_j[i].z);
    }
}

/*
 * @fn energy()
 *
 * @brief
 *  Kinetic energy calculation.
 *
 *  The Potential energy is calculated in the init_acc_jrk()
 *    and in the update_acc_jrk() functions.
 *
 */
float energy()
{
    double ekin_tmp = 0;
    double epot_tmp = 0;
    double r2, v2;
    double rx, ry, rz;
    double vx, vy, vz;
    int i, j;

    epot = 0;
    ekin = 0;

//    #pragma omp parallel for private(epot_tmp,j,rx,ry,rz,r2,vx,vy,vz,v2,ekin_tmp)
    for (i = 0; i < n; i++)
    {
        epot_tmp = 0;
        for (j = i+1; j < n; j++)
        {
            rx = h_r[j].x - h_r[i].x;
            ry = h_r[j].y - h_r[i].y;
            rz = h_r[j].z - h_r[i].z;

            r2 = rx*rx + ry*ry + rz*rz;

            epot_tmp -= (h_m[i] * h_m[j]) / sqrt(r2);
        }

        vx = h_v[i].x * h_v[i].x;
        vy = h_v[i].y * h_v[i].y;
        vz = h_v[i].z * h_v[i].z;

        v2 = vx + vy + vz;

        ekin_tmp = 0.5 * h_m[i] * v2;

        ekin += ekin_tmp;
        epot += epot_tmp;
    }
    return epot + ekin;
}

void get_energy_log(double ITIME, int iterations, int nsteps)
{
    energy_end = energy();
    double relative_error = abs((energy_end-energy_ini)/energy_ini);
    energy_tmp += relative_error;

    if((int)ITIME == 1)
        printf("ITIME  ITER NSTEPS/ITER    ETA     TIME            ENERGY     REL_ENER     CUM_ENER\n");

    printf("%4d %6d %.6f %f %.6f %.15e %.15e %.15e\n",
            (int)ITIME,
            iterations,
            nsteps/(float)iterations,
            ETA_N,
            omp_get_wtime() - ini_time,
            energy_end,
            relative_error,
            energy_tmp/(int)ITIME);
}

/*
 * @fn save_old()
 *
 * @brief
 *  Save the previous particle position, velocity,
 *    acceleration and jrk.
 */
void save_old(int total)
{
    for (int k = 0; k < total; k++)
    {
        int i = h_move[k];
        h_old_a[i].x = h_a[i].x;
        h_old_a[i].y = h_a[i].y;
        h_old_a[i].z = h_a[i].z;

        h_old_j[i].x = h_j[i].x;
        h_old_j[i].y = h_j[i].y;
        h_old_j[i].z = h_j[i].z;
    }
}

/*
 * @fn predicted_pos_vel()
 *
 * @brief
 *  Perform the calculation of the predicted position
 *    and velocity, following the Hermite integrator scheme
 *    formulas.
 */
void
predicted_pos_vel(double ITIME)
{
    //#pragma omp parallel for
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double dt = ITIME - h_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        h_p_r[i].x = (dt3/6 * h_j[i].x) + (dt2/2 * h_a[i].x) + (dt * h_v[i].x) + h_r[i].x;
        h_p_r[i].y = (dt3/6 * h_j[i].y) + (dt2/2 * h_a[i].y) + (dt * h_v[i].y) + h_r[i].y;
        h_p_r[i].z = (dt3/6 * h_j[i].z) + (dt2/2 * h_a[i].z) + (dt * h_v[i].z) + h_r[i].z;

        h_p_v[i].x = (dt2/2 * h_j[i].x) + (dt * h_a[i].x) + h_v[i].x;
        h_p_v[i].y = (dt2/2 * h_j[i].y) + (dt * h_a[i].y) + h_v[i].y;
        h_p_v[i].z = (dt2/2 * h_j[i].z) + (dt * h_a[i].z) + h_v[i].z;
        //h_p_r[i].x += dt3 * h_j[i].x + dt2 * h_a[i].x + dt * h_v[i].x + h_r[i].x;
        //h_p_r[i].y += dt3 * h_j[i].y + dt2 * h_a[i].y + dt * h_v[i].y + h_r[i].y;
        //h_p_r[i].z += dt3 * h_j[i].z + dt2 * h_a[i].z + dt * h_v[i].z + h_r[i].z;

        //h_p_v[i].x += dt2 * h_j[i].x + dt * h_a[i].x + h_v[i].x;
        //h_p_v[i].y += dt2 * h_j[i].y + dt * h_a[i].y + h_v[i].y;
        //h_p_v[i].z += dt2 * h_j[i].z + dt * h_a[i].z + h_v[i].z;
    }
}

/*
 * @fn predicted_pos_vel_kepler()
 *
 * @brief
 *  Perform the calculation of the predicted position
 *    and velocity using the Kepler's equation.
 *
 */
void predicted_pos_vel_kepler(float ITIME, int total)
{

    for (int i = 0; i < total; i++)
    {
        int k = h_move[i];
        float dt = ITIME - h_t[i];

        double rx = h_r[k].x - h_r[0].x;
        double ry = h_r[k].y - h_r[0].y;
        double rz = h_r[k].z - h_r[0].z;

        double vx = h_v[k].x - h_v[0].x;
        double vy = h_v[k].y - h_v[0].y;
        double vz = h_v[k].z - h_v[0].z;

        for (float j = 0; j < KEPLER_ITE; j+=dt)
        {
            kepler_prediction(&rx, &ry, &rz, &vx, &vy, &vz, h_m[k], h_m[0], j, k);
        }
        h_p_r[k].x = rx;
        h_p_r[k].y = ry;
        h_p_r[k].z = rz;

        h_p_v[k].x = vx;
        h_p_v[k].y = vy;
        h_p_v[k].z = vz;
    }
}

/*
 * @fn correction_pos_vel()
 *
 * @brief
 *  Perform the calculation of the acceleration derivates
 *    to apply the correction terms to the position and velocity.
 *
 */
void correction_pos_vel(double ITIME, int total)
{

    for (int k = 0; k < total; k++)
    {
        int i = h_move[k];

        double dt1 = h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        double ax0_2 = (-6 * (h_old_a[i].x - h_a[i].x ) - dt1 * (4 * h_old_j[i].x + 2 * h_j[i].x) ) / dt2;
        double ay0_2 = (-6 * (h_old_a[i].y - h_a[i].y ) - dt1 * (4 * h_old_j[i].y + 2 * h_j[i].y) ) / dt2;
        double az0_2 = (-6 * (h_old_a[i].z - h_a[i].z ) - dt1 * (4 * h_old_j[i].z + 2 * h_j[i].z) ) / dt2;

        // Acceleration 3rd derivate
        double ax0_3 = (12 * (h_old_a[i].x - h_a[i].x ) + 6 * dt1 * (h_old_j[i].x + h_j[i].x) ) / dt3;
        double ay0_3 = (12 * (h_old_a[i].y - h_a[i].y ) + 6 * dt1 * (h_old_j[i].y + h_j[i].y) ) / dt3;
        double az0_3 = (12 * (h_old_a[i].z - h_a[i].z ) + 6 * dt1 * (h_old_j[i].z + h_j[i].z) ) / dt3;

        // Correcting position
        h_r[i].x = h_p_r[i].x + (dt4/24)*ax0_2 + (dt5/120)*ax0_3;
        h_r[i].y = h_p_r[i].y + (dt4/24)*ay0_2 + (dt5/120)*ay0_3;
        h_r[i].z = h_p_r[i].z + (dt4/24)*az0_2 + (dt5/120)*az0_3;

        // Correcting velocity
        h_v[i].x = h_p_v[i].x + (dt3/6)*ax0_2 + (dt4/24)*ax0_3;
        h_v[i].y = h_p_v[i].y + (dt3/6)*ay0_2 + (dt4/24)*ay0_3;
        h_v[i].z = h_p_v[i].z + (dt3/6)*az0_2 + (dt4/24)*az0_3;

        // Timestep update

        // Calculating a_{1,i}^{(2)} = a_{0,i}^{(2)} + dt * a_{0,i}^{(3)}
        double ax1_2 = ax0_2 + dt1 * ax0_3;
        double ay1_2 = ay0_2 + dt1 * ay0_3;
        double az1_2 = az0_2 + dt1 * az0_3;

        // |a_{1,i}|
        double abs_a1 = get_magnitude(h_a[i].x, h_a[i].y, h_a[i].z);
        // |j_{1,i}|
        double abs_j1 = get_magnitude(h_j[i].x, h_j[i].y, h_j[i].z);
        // |j_{1,i}|^{2}
        double abs_j12  = abs_j1 * abs_j1;
        // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
        double abs_a1_3 = get_magnitude(ax0_3, ay0_3, az0_3);
        // |a_{1,i}^{(2)}|
        double abs_a1_2 = get_magnitude(ax1_2, ay1_2, az1_2);
        // |a_{1,i}^{(2)}|^{2}
        double abs_a1_22  = abs_a1_2 * abs_a1_2;

        double tmp_dt = sqrt(ETA_N * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));

        h_t[i] = ITIME;
        tmp_dt = normalize_dt(tmp_dt, h_dt[i], h_t[i], i);
        h_dt[i] = tmp_dt;
    }
}

double normalize_dt(double new_dt, double old_dt, double t, int i)
{
    if (new_dt <= old_dt/8)
    {
        new_dt = D_TIME_MIN;
    }
    else if ( old_dt/8 < new_dt && new_dt <= old_dt/4)
    {
        new_dt = old_dt / 8;
    }
    else if ( old_dt/4 < new_dt && new_dt <= old_dt/2)
    {
        new_dt = old_dt / 4;
    }
    else if ( old_dt/2 < new_dt && new_dt <= old_dt)
    {
        new_dt = old_dt / 2;
    }
    else if ( old_dt < new_dt && new_dt <= old_dt * 2)
    {
        new_dt = old_dt;
    }
    else if (2 * old_dt < new_dt)
    {
        float val = t/(2 *old_dt);
        if(std::ceil(val) == val)
        {
            new_dt = 2 * old_dt;
        }
        else
        {
            new_dt = old_dt;
        }
    }
    else
    {
        std::cerr << "Nothing to do!" << std::endl;
        getchar();
    }

    if (new_dt < D_TIME_MIN)
    {
        new_dt = D_TIME_MIN;
    }
    else if (new_dt > D_TIME_MAX)
    {
        new_dt = D_TIME_MAX;
    }

    return new_dt;
}

//    Point p = get_center_of_density();
//    for (int i = 0; i < n; i++) {
//        double rx = h_r[i].x - p.x;
//        double ry = h_r[i].y - p.y;
//        double rz = h_r[i].z - p.z;
//        double d  = get_magnitude(rx, ry, rz);
//        printf("%d %.6f\n",i, d);
//    }

//    float t_rh = get_relaxation_time();
//    float t_cr = get_crossing_time();
//
//    std::cout << "T_rh : " << t_rh << std::endl;
//    std::cout << "T_cr : " << t_cr << std::endl;
//    std::cout << "T_cc : " << 17 * t_rh << std::endl;
//
