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
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] next_itime()\n");
    #endif
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
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] find_particles_to_move()\n");
    #endif
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        h_move[i] = -1;
        if (h_t[i] + h_dt[i] == ITIME)
        {
            #ifdef USE_KEPLER
            if(i != 0)
            {
                h_move[j] = i;
                j++;
            }
            #else
                h_move[j] = i;
                j++;
            #endif
//            std::cout << i << " ";
        }
    }
 //   std::cout << std::endl;
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
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] init_dt()\n");
    #endif
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    double tmp_dt;
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        #ifdef USE_KEPLER
        if(i == 0)
        {
            h_dt[0] = 0.0;
            continue;
        }
        #endif
        double a2 = get_magnitude(h_a[i].x, h_a[i].y, h_a[i].z);
        double j2 = get_magnitude(h_a1[i].x, h_a1[i].y, h_a1[i].z);
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

void init_dt2(double *ATIME)
{
    // Get Snap and Crackle
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] init_dt2()\n");
    #endif
    int i = 0;
    int j = 0;
    for (i = INIT_PARTICLE; i < n; i++) {
        for (j = INIT_PARTICLE; i < n; i++) {
            double rx = h_r[j].x - h_r[i].x;
            double ry = h_r[j].y - h_r[i].y;
            double rz = h_r[j].z - h_r[i].z;

            double vx = h_v[j].x - h_v[i].x;
            double vy = h_v[j].y - h_v[i].y;
            double vz = h_v[j].z - h_v[i].z;

            double r2 = rx*rx + ry*ry + rz*rz + softening*softening;
            double rinv = 1/sqrt(r2);
            double r2inv = rinv  * rinv;
            double r3inv = r2inv * rinv;
            double mr3inv = r3inv * h_m[j];

            double rv = rx*vx + ry*vy + rz*vz;
            double ra = rx*h_a[i].x + ry*h_a[i].y + rz*h_a[i].z;
            double v2 = vx*vx + vy*vy + vz*vz;
            double va = h_a[i].x*vx + h_a[i].y*vy + h_a[i].z*vz;
            double rj = rx*h_a1[i].x + ry*h_a1[i].y + rz*h_a1[i].z;

            double alpha = rv * r2inv;
            double beta = alpha*alpha + r2inv * (v2 + ra);
            double gamma = 3 * va + rj * r2inv + alpha * (3*beta - 4*alpha*alpha);

            h_a2[i].x += mr3inv*h_a[i].x - 6*alpha*h_a1[i].x - 3*beta*h_a[i].x;
            h_a2[i].y += mr3inv*h_a[i].y - 6*alpha*h_a1[i].y - 3*beta*h_a[i].y;
            h_a2[i].z += mr3inv*h_a[i].z - 6*alpha*h_a1[i].z - 3*beta*h_a[i].z;

            h_a3[i].x += mr3inv*h_a1[i].x - 9*alpha*h_a2[i].x - 9*beta*h_a1[i].x -3*gamma*h_a[i].x;
            h_a3[i].y += mr3inv*h_a1[i].y - 9*alpha*h_a2[i].y - 9*beta*h_a1[i].y -3*gamma*h_a[i].y;
            h_a3[i].z += mr3inv*h_a1[i].z - 9*alpha*h_a2[i].z - 9*beta*h_a1[i].z -3*gamma*h_a[i].z;
        }
    }
    // Get Timesteps
    for (i = 0; i < n; i++) {
        #ifdef USE_KEPLER
        if(i == 0)
        {
            h_dt[0] = 0.0;
            continue;
        }
        #endif
        double m_a = get_magnitude(h_a[i].x, h_a[i].y, h_a[i].z);
        double m_a1 = get_magnitude(h_a1[i].x, h_a1[i].y, h_a1[i].z);
        double m_a2 = get_magnitude(h_a2[i].x, h_a2[i].y, h_a2[i].z);
        double m_a3 = get_magnitude(h_a3[i].x, h_a3[i].y, h_a3[i].z);
        double tmp_dt = sqrt(ETA_S * (m_a * m_a2 + m_a1*m_a1)/(m_a1*m_a3+m_a2*m_a2));
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

    double r2 = rx*rx + ry*ry + rz*rz + softening*softening;
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

    h_a1[i].x += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    h_a1[i].y += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    h_a1[i].z += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void init_acc_jrk()
{
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] init_acc_jrk()\n");
    #endif
    int j;
    #pragma omp parallel for private(j)
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
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] update_acc_jrk()\n");
    #endif
    int i, j;
    for (int k = 0; k < total; k++)
    {
        i = h_move[k];
        h_a[i].x = 0.0;
        h_a[i].y = 0.0;
        h_a[i].z = 0.0;
        h_a1[i].x = 0.0;
        h_a1[i].y = 0.0;
        h_a1[i].z = 0.0;

        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i,j);
        }
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
double energy()
{
    double ekin_tmp;
    double epot_tmp;
    int i, j;

    epot = 0.0;
    ekin = 0.0;

    for (i = 0; i < n; i++)
    {
        epot_tmp = 0.0;
        for (j = i+1; j < n; j++)
        {
            double rx = h_r[j].x - h_r[i].x;
            double ry = h_r[j].y - h_r[i].y;
            double rz = h_r[j].z - h_r[i].z;
            double r2 = rx*rx + ry*ry + rz*rz;
            epot_tmp -= (h_m[i] * h_m[j]) / sqrt(r2);
        }

        double vx = h_v[i].x * h_v[i].x;
        double vy = h_v[i].y * h_v[i].y;
        double vz = h_v[i].z * h_v[i].z;
        double v2 = vx + vy + vz;
        ekin_tmp = 0.5 * h_m[i] * v2;

        ekin += ekin_tmp;
        epot += epot_tmp;
    }
    return epot + ekin;
}

void get_energy_log(double ITIME, int iterations, int nsteps, FILE *out)
{
    energy_end = energy();
    double relative_error   = abs((energy_end-energy_tmp)/energy_ini);
    double cumulative_error = abs((energy_end-energy_ini)/energy_ini);
    energy_tmp = energy_end;

    printf("% 6d % 10d % 10d % .6f % .6f % .6f % .15e % .15e % .15e 00\n",
    //fprintf(out, "% 6d % 10d % 10d % .6f % .6f % .6f % .15e % .15e % .15e 00\n",
            (int)ITIME,
            iterations,
            nsteps,
            eta,
            softening,
            (float)clock()/CLOCKS_PER_SEC - ini_time,
            energy_end,
            relative_error,
            cumulative_error);
    fflush(out);
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
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] save_old()\n");
    #endif
    for (int k = INIT_PARTICLE; k < total; k++)
    {
        int i = h_move[k];
        h_old_a[i].x = h_a[i].x;
        h_old_a[i].y = h_a[i].y;
        h_old_a[i].z = h_a[i].z;

        h_old_a1[i].x = h_a1[i].x;
        h_old_a1[i].y = h_a1[i].y;
        h_old_a1[i].z = h_a1[i].z;
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
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] predicted_pos_vel()\n");
    #endif
    //#pragma omp parallel for
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double dt = ITIME - h_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        #ifdef USE_KEPLER
        if (h_move[i] == -1)
        {
            h_p_r[i].x = (dt3/6 * h_a1[i].x) + (dt2/2 * h_a[i].x) + (dt * h_v[i].x) + h_r[i].x;
            h_p_r[i].y = (dt3/6 * h_a1[i].y) + (dt2/2 * h_a[i].y) + (dt * h_v[i].y) + h_r[i].y;
            h_p_r[i].z = (dt3/6 * h_a1[i].z) + (dt2/2 * h_a[i].z) + (dt * h_v[i].z) + h_r[i].z;

            h_p_v[i].x = (dt2/2 * h_a1[i].x) + (dt * h_a[i].x) + h_v[i].x;
            h_p_v[i].y = (dt2/2 * h_a1[i].y) + (dt * h_a[i].y) + h_v[i].y;
            h_p_v[i].z = (dt2/2 * h_a1[i].z) + (dt * h_a[i].z) + h_v[i].z;
        }
        else
        {
            h_p_r[i].x += (dt3/6 * h_a1[i].x) + (dt2/2 * h_a[i].x);
            h_p_r[i].y += (dt3/6 * h_a1[i].y) + (dt2/2 * h_a[i].y);
            h_p_r[i].z += (dt3/6 * h_a1[i].z) + (dt2/2 * h_a[i].z);

            h_p_v[i].x += (dt2/6 * h_a1[i].x) + (dt/2 * h_a[i].x);
            h_p_v[i].y += (dt2/6 * h_a1[i].y) + (dt/2 * h_a[i].y);
            h_p_v[i].z += (dt2/6 * h_a1[i].z) + (dt/2 * h_a[i].z);

        }
        #else
        h_p_r[i].x = (dt3/6 * h_a1[i].x) + (dt2/2 * h_a[i].x) + (dt * h_v[i].x) + h_r[i].x;
        h_p_r[i].y = (dt3/6 * h_a1[i].y) + (dt2/2 * h_a[i].y) + (dt * h_v[i].y) + h_r[i].y;
        h_p_r[i].z = (dt3/6 * h_a1[i].z) + (dt2/2 * h_a[i].z) + (dt * h_v[i].z) + h_r[i].z;

        h_p_v[i].x = (dt2/2 * h_a1[i].x) + (dt * h_a[i].x) + h_v[i].x;
        h_p_v[i].y = (dt2/2 * h_a1[i].y) + (dt * h_a[i].y) + h_v[i].y;
        h_p_v[i].z = (dt2/2 * h_a1[i].z) + (dt * h_a[i].z) + h_v[i].z;
        #endif

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
void predicted_pos_vel_kepler(double ITIME, int total)
{
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] predicted_pos_vel_kepler()\n");
    #endif
    for (int i = 0; i < total; i++)
    {
        int k = h_move[i];
        double dt = ITIME - h_t[k];
        //double time = 0.0;


        double rx = h_r[k].x - h_r[0].x;
        double ry = h_r[k].y - h_r[0].y;
        double rz = h_r[k].z - h_r[0].z;

        double vx = h_v[k].x - h_v[0].x;
        double vy = h_v[k].y - h_v[0].y;
        double vz = h_v[k].z - h_v[0].z;
        #ifdef DEBUG_KEPLER
        printf("Particle %d\n", k);
        printf("[Old position] %.15f %.15f %.15f\n", h_r[k].x, h_r[k].y, h_r[k].z);
        printf("[Old velocity] %.15f %.15f %.15f\n", h_v[k].x, h_v[k].y, h_v[k].z);
        #endif

        //for (time = h_t[k]; time < ITIME; time+=dt)
        //{
        //    kepler_prediction(&rx, &ry, &rz, &vx, &vy, &vz, dt, k);
        //}
        kepler_prediction(&rx, &ry, &rz, &vx, &vy, &vz, dt, k);

        h_p_r[k].x = rx;
        h_p_r[k].y = ry;
        h_p_r[k].z = rz;

        h_p_v[k].x = vx;
        h_p_v[k].y = vy;
        h_p_v[k].z = vz;
        #ifdef DEBUG_KEPLER
        printf("[New position] %.15f %.15f %.15f\n", h_p_r[k].x, h_p_r[k].y, h_p_r[k].z);
        printf("[New velocity] %.15f %.15f %.15f\n", h_p_v[k].x, h_p_v[k].y, h_p_v[k].z);
        printf("End particle %d\n", k);
        getchar();
        #endif
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
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] correction_pos_vel()\n");
    #endif

    for (int k = 0; k < total; k++)
    {
        int i = h_move[k];

        double dt1 = h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        #ifdef USE_KEPLER
        h_a2[i].x += (-6 * (h_old_a[i].x - h_a[i].x ) - dt1 * (4 * h_old_a1[i].x + 2 * h_a1[i].x) ) / dt2;
        h_a2[i].y += (-6 * (h_old_a[i].y - h_a[i].y ) - dt1 * (4 * h_old_a1[i].y + 2 * h_a1[i].y) ) / dt2;
        h_a2[i].z += (-6 * (h_old_a[i].z - h_a[i].z ) - dt1 * (4 * h_old_a1[i].z + 2 * h_a1[i].z) ) / dt2;
        #else
        h_a2[i].x = (-6 * (h_old_a[i].x - h_a[i].x ) - dt1 * (4 * h_old_a1[i].x + 2 * h_a1[i].x) ) / dt2;
        h_a2[i].y = (-6 * (h_old_a[i].y - h_a[i].y ) - dt1 * (4 * h_old_a1[i].y + 2 * h_a1[i].y) ) / dt2;
        h_a2[i].z = (-6 * (h_old_a[i].z - h_a[i].z ) - dt1 * (4 * h_old_a1[i].z + 2 * h_a1[i].z) ) / dt2;
        #endif

        // Acceleration 3rd derivate
        #ifdef USE_KEPLER
        h_a3[i].x += (12 * (h_old_a[i].x - h_a[i].x ) + 6 * dt1 * (h_old_a1[i].x + h_a1[i].x) ) / dt3;
        h_a3[i].y += (12 * (h_old_a[i].y - h_a[i].y ) + 6 * dt1 * (h_old_a1[i].y + h_a1[i].y) ) / dt3;
        h_a3[i].z += (12 * (h_old_a[i].z - h_a[i].z ) + 6 * dt1 * (h_old_a1[i].z + h_a1[i].z) ) / dt3;
        #else
        h_a3[i].x = (12 * (h_old_a[i].x - h_a[i].x ) + 6 * dt1 * (h_old_a1[i].x + h_a1[i].x) ) / dt3;
        h_a3[i].y = (12 * (h_old_a[i].y - h_a[i].y ) + 6 * dt1 * (h_old_a1[i].y + h_a1[i].y) ) / dt3;
        h_a3[i].z = (12 * (h_old_a[i].z - h_a[i].z ) + 6 * dt1 * (h_old_a1[i].z + h_a1[i].z) ) / dt3;
        #endif

        // Correcting position
        h_r[i].x = h_p_r[i].x + (dt4/24)*h_a2[i].x + (dt5/120)*h_a3[i].x;
        h_r[i].y = h_p_r[i].y + (dt4/24)*h_a2[i].y + (dt5/120)*h_a3[i].y;
        h_r[i].z = h_p_r[i].z + (dt4/24)*h_a2[i].z + (dt5/120)*h_a3[i].z;

        // Correcting velocity
        h_v[i].x = h_p_v[i].x + (dt3/6)*h_a2[i].x + (dt4/24)*h_a3[i].x;
        h_v[i].y = h_p_v[i].y + (dt3/6)*h_a2[i].y + (dt4/24)*h_a3[i].y;
        h_v[i].z = h_p_v[i].z + (dt3/6)*h_a2[i].z + (dt4/24)*h_a3[i].z;

        h_t[i] = ITIME;
        double normal_dt  = get_timestep_normal(i);
        #ifdef USE_KEPLER
        double central_dt = get_timestep_central(i);
        if(normal_dt > central_dt)
        {
            central_dt = normalize_dt(central_dt, h_dt[i], h_t[i], i);
            h_dt[i] = central_dt;
        }
        else
        {
            normal_dt = normalize_dt(normal_dt, h_dt[i], h_t[i], i);
            h_dt[i] = normal_dt;
        }
        #else
            normal_dt = normalize_dt(normal_dt, h_dt[i], h_t[i], i);
            h_dt[i] = normal_dt;
        #endif

    }
}

double get_timestep_normal(int i)
{
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] get_timestep_normal()\n");
    #endif
    // Calculating a_{1,i}^{(2)} = a_{0,i}^{(2)} + dt * a_{0,i}^{(3)}
    double ax1_2 = h_a2[i].x + h_dt[i] * h_a3[i].x;
    double ay1_2 = h_a2[i].y + h_dt[i] * h_a3[i].y;
    double az1_2 = h_a2[i].z + h_dt[i] * h_a3[i].z;

    // |a_{1,i}|
    double abs_a1 = get_magnitude(h_a[i].x, h_a[i].y, h_a[i].z);
    // |j_{1,i}|
    double abs_j1 = get_magnitude(h_a1[i].x, h_a1[i].y, h_a1[i].z);
    // |j_{1,i}|^{2}
    double abs_j12  = abs_j1 * abs_j1;
    // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
    double abs_a1_3 = get_magnitude(h_a3[i].x, h_a3[i].y, h_a3[i].z);
    // |a_{1,i}^{(2)}|
    double abs_a1_2 = get_magnitude(ax1_2, ay1_2, az1_2);
    // |a_{1,i}^{(2)}|^{2}
    double abs_a1_22  = abs_a1_2 * abs_a1_2;

    double normal_dt = sqrt(eta * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));


    return normal_dt;
}

double get_timestep_central(int i)
{
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] get_timestep_central()\n");
    #endif
    double r = get_magnitude(h_r[i].x, h_r[i].y, h_r[i].z);
    double r3 = r*r*r;
    double central_dt = ((2.0 * M_PI )/OSTEPS * sqrt(r3/(G * h_m[0])));

    return central_dt;
}

double normalize_dt(double new_dt, double old_dt, double t, int i)
{
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] normalize_dt()\n");
    #endif
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
    //else
    //{
    //    //fprintf(stderr, "gravidy: Undefined treatment for the time-step of (%d)",i);
    //    new_dt = old_dt;
    //}

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
