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
            //printf("%d ",i);
            //printf("pre ");
            //print_particle(i);
        }
    }
    //printf("\n");
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
        double a2 = get_magnitude(h_f[i].a[0], h_f[i].a[1],h_f[i].a[2]);
        double j2 = get_magnitude(h_f[i].a1[0], h_f[i].a1[1],h_f[i].a1[2]);
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
    double rx = h_p[j].r[0] - h_p[i].r[0];
    double ry = h_p[j].r[1] - h_p[i].r[1];
    double rz = h_p[j].r[2] - h_p[i].r[2];

    double vx = h_p[j].v[0] - h_p[i].v[0];
    double vy = h_p[j].v[1] - h_p[i].v[1];
    double vz = h_p[j].v[2] - h_p[i].v[2];

    double r2     = rx*rx + ry*ry + rz*rz + e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * h_m[j];
    double mr5inv = r5inv * h_m[j];

    double rv = rx*vx + ry*vy + rz*vz;

    h_f[i].a[0] += (rx * mr3inv);
    h_f[i].a[1] += (ry * mr3inv);
    h_f[i].a[2] += (rz * mr3inv);

    h_f[i].a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    h_f[i].a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    h_f[i].a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void init_acc_jrk()
{
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] init_acc_jrk()\n");
    #endif
    //#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d %d\n",i, j);
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
    for (int k = 0; k < total; k++)
    {
        int i = h_move[k];
        h_f[i].a[0]  = 0.0;
        h_f[i].a[1]  = 0.0;
        h_f[i].a[2]  = 0.0;
        h_f[i].a1[0] = 0.0;
        h_f[i].a1[1] = 0.0;
        h_f[i].a1[2] = 0.0;

        for (int j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i,j);
        }
//        printf("Updating %d %f\t%f\t%f\t%f\t%f\t%f\n",
//                i, h_f[k].a[0],  h_f[k].a[1], h_f[k].a[2],
//                   h_f[k].a1[0], h_f[k].a[1], h_f[k].a[2]);
//
//        getchar();
    }
    gtime.update_end += omp_get_wtime() - gtime.update_ini;
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
    epot = 0.0;
    ekin = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double epot_tmp = 0.0;
        for (int j = i+1; j < n; j++)
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

        double ekin_tmp = 0.5 * h_m[i] * v2;

        #pragma omp atomic
        ekin += ekin_tmp;
        #pragma omp atomic
        epot += epot_tmp;
    }
    return epot + ekin;
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
        h_old_a[i].x = h_f[i].a[0];
        h_old_a[i].y = h_f[i].a[1];
        h_old_a[i].z = h_f[i].a[2];

        h_old_a1[i].x = h_f[i].a1[0];
        h_old_a1[i].y = h_f[i].a1[1];
        h_old_a1[i].z = h_f[i].a1[2];
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
void predicted_pos_vel(double ITIME)
{
    #ifdef DEBUG_HERMITE
    printf("[DEBUG] predicted_pos_vel()\n");
    #endif
    gtime.prediction_ini = omp_get_wtime();
    //#pragma omp parallel for
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double dt  = ITIME - h_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        #ifdef USE_KEPLER
        if (h_move[i] == -1)
        {

            h_p[i].r[0] = (dt3/6 * h_f[i].a1[0]) + (dt2/2 * h_f[i].a[0]) + (dt * h_v[i].x) + h_r[i].x;
            h_p[i].r[1] = (dt3/6 * h_f[i].a1[1]) + (dt2/2 * h_f[i].a[1]) + (dt * h_v[i].y) + h_r[i].y;
            h_p[i].r[2] = (dt3/6 * h_f[i].a1[2]) + (dt2/2 * h_f[i].a[2]) + (dt * h_v[i].z) + h_r[i].z;

            h_p[i].v[0] = (dt2/2 * h_f[i].a1[0]) + (dt * h_f[i].a[0]) + h_v[i].x;
            h_p[i].v[1] = (dt2/2 * h_f[i].a1[1]) + (dt * h_f[i].a[1]) + h_v[i].y;
            h_p[i].v[2] = (dt2/2 * h_f[i].a1[2]) + (dt * h_f[i].a[2]) + h_v[i].z;
        }
        else
        {
            h_p[i].r[0] += (dt3/6 * h_f[i].a1[0]) + (dt2/2 * h_f[i].a[0]);
            h_p[i].r[1] += (dt3/6 * h_f[i].a1[1]) + (dt2/2 * h_f[i].a[1]);
            h_p[i].r[2] += (dt3/6 * h_f[i].a1[2]) + (dt2/2 * h_f[i].a[2]);

            h_p[i].v[0] += (dt2/2 * h_f[i].a1[0]) + (dt * h_f[i].a[0]);
            h_p[i].v[1] += (dt2/2 * h_f[i].a1[1]) + (dt * h_f[i].a[1]);
            h_p[i].v[2] += (dt2/2 * h_f[i].a1[2]) + (dt * h_f[i].a[2]);

        }
        #else
        h_p[i].r[0] = (dt3/6 * h_f[i].a1[0]) + (dt2/2 * h_f[i].a[0]) + (dt * h_v[i].x) + h_r[i].x;
        h_p[i].r[1] = (dt3/6 * h_f[i].a1[1]) + (dt2/2 * h_f[i].a[1]) + (dt * h_v[i].y) + h_r[i].y;
        h_p[i].r[2] = (dt3/6 * h_f[i].a1[2]) + (dt2/2 * h_f[i].a[2]) + (dt * h_v[i].z) + h_r[i].z;

        h_p[i].v[0] = (dt2/2 * h_f[i].a1[0]) + (dt * h_f[i].a[0]) + h_v[i].x;
        h_p[i].v[1] = (dt2/2 * h_f[i].a1[1]) + (dt * h_f[i].a[1]) + h_v[i].y;
        h_p[i].v[2] = (dt2/2 * h_f[i].a1[2]) + (dt * h_f[i].a[2]) + h_v[i].z;
        #endif

        //printf("pred_normal %d\n",i);
        //print_particle(i);
        //getchar();

    }
    gtime.prediction_end += omp_get_wtime() - gtime.prediction_ini;
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
        printf("[Old acceleration] %.15f %.15f %.15f\n", h_f[k].a[0], h_f[k].a[1], h_f[k].a[2]);
        printf("[Old jerk] %.15f %.15f %.15f\n", h_f[k].a1[0], h_f[k].a1[1], h_f[k].a1[2]);
        #endif

        //for (time = h_t[k]; time < ITIME; time+=dt)
        //{
        //    kepler_prediction(&rx, &ry, &rz, &vx, &vy, &vz, dt, k);
        //}
        kepler_prediction(&rx, &ry, &rz, &vx, &vy, &vz, dt, k);

        h_p[k].r[0] = rx;
        h_p[k].r[1] = ry;
        h_p[k].r[2] = rz;

        h_p[k].v[0] = vx;
        h_p[k].v[1] = vy;
        h_p[k].v[2] = vz;
        #ifdef DEBUG_KEPLER
        printf("[New position] %.15f %.15f %.15f\n", h_p[k].r[0], h_p[k].r[1], h_p[k].r[2]);
        printf("[New velocity] %.15f %.15f %.15f\n", h_p[k].v[0], h_p[k].v[1], h_p[k].v[2]);
        printf("[New acceleration] %.15f %.15f %.15f\n", h_f[k].a[0], h_f[k].a[1], h_f[k].a[2]);
        printf("[New jerk] %.15f %.15f %.15f\n", h_f[k].a1[0], h_f[k].a1[1], h_f[k].a1[2]);
        printf("End particle %d\n", k);
        getchar();
        #endif

        //printf("pred_kepler ");
        //print_particle(k);
        //getchar();
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
    gtime.correction_ini = omp_get_wtime();
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
        h_a2[i].x += (-6 * (h_old_a[i].x - h_f[i].a[0] ) - dt1 * (4 * h_old_a1[i].x + 2 * h_f[i].a1[0]) ) / dt2;
        h_a2[i].y += (-6 * (h_old_a[i].y - h_f[i].a[1] ) - dt1 * (4 * h_old_a1[i].y + 2 * h_f[i].a1[1]) ) / dt2;
        h_a2[i].z += (-6 * (h_old_a[i].z - h_f[i].a[2] ) - dt1 * (4 * h_old_a1[i].z + 2 * h_f[i].a1[2]) ) / dt2;
        #else
        h_a2[i].x = (-6 * (h_old_a[i].x - h_f[i].a[0] ) - dt1 * (4 * h_old_a1[i].x + 2 * h_f[i].a1[0]) ) / dt2;
        h_a2[i].y = (-6 * (h_old_a[i].y - h_f[i].a[1] ) - dt1 * (4 * h_old_a1[i].y + 2 * h_f[i].a1[1]) ) / dt2;
        h_a2[i].z = (-6 * (h_old_a[i].z - h_f[i].a[2] ) - dt1 * (4 * h_old_a1[i].z + 2 * h_f[i].a1[2]) ) / dt2;
        #endif

        // Acceleration 3rd derivate
        #ifdef USE_KEPLER
        h_a3[i].x += (12 * (h_old_a[i].x - h_f[i].a[0] ) + 6 * dt1 * (h_old_a1[i].x + h_f[i].a1[0]) ) / dt3;
        h_a3[i].y += (12 * (h_old_a[i].y - h_f[i].a[1] ) + 6 * dt1 * (h_old_a1[i].y + h_f[i].a1[1]) ) / dt3;
        h_a3[i].z += (12 * (h_old_a[i].z - h_f[i].a[2] ) + 6 * dt1 * (h_old_a1[i].z + h_f[i].a1[2]) ) / dt3;
        #else
        h_a3[i].x = (12 * (h_old_a[i].x - h_f[i].a[0] ) + 6 * dt1 * (h_old_a1[i].x + h_f[i].a1[0]) ) / dt3;
        h_a3[i].y = (12 * (h_old_a[i].y - h_f[i].a[1] ) + 6 * dt1 * (h_old_a1[i].y + h_f[i].a1[1]) ) / dt3;
        h_a3[i].z = (12 * (h_old_a[i].z - h_f[i].a[2] ) + 6 * dt1 * (h_old_a1[i].z + h_f[i].a1[2]) ) / dt3;
        #endif

        // Correcting position
        h_r[i].x = h_p[i].r[0] + (dt4/24)*h_a2[i].x + (dt5/120)*h_a3[i].x;
        h_r[i].y = h_p[i].r[1] + (dt4/24)*h_a2[i].y + (dt5/120)*h_a3[i].y;
        h_r[i].z = h_p[i].r[2] + (dt4/24)*h_a2[i].z + (dt5/120)*h_a3[i].z;

        // Correcting velocity
        h_v[i].x = h_p[i].v[0] + (dt3/6)*h_a2[i].x + (dt4/24)*h_a3[i].x;
        h_v[i].y = h_p[i].v[1] + (dt3/6)*h_a2[i].y + (dt4/24)*h_a3[i].y;
        h_v[i].z = h_p[i].v[2] + (dt3/6)*h_a2[i].z + (dt4/24)*h_a3[i].z;


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

        //printf("correction ");
        //print_particle(i);
        //getchar();

    }


    gtime.correction_end += omp_get_wtime() - gtime.correction_ini;
}

/*
 * @fn get_timestep_normal()
 *
 * @brief
 *  Obtain the time step of an i particle
 *
 */
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
    double abs_a1 = get_magnitude(h_f[i].a[0], h_f[i].a[1], h_f[i].a[2]);
    // |j_{1,i}|
    double abs_j1 = get_magnitude(h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2]);
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
    double r = get_magnitude(h_r[i].x-h_r[0].x, h_r[i].y-h_r[0].y, h_r[i].z-h_r[i].z);
    double r3 = r*r*r;
    double central_dt = (((2.0 * M_PI )/OSTEPS) * sqrt(r3/(G * h_m[0])));

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
    else
    {
        fprintf(stderr, "gravidy: Undefined treatment for the time-step of (%d)",i);
        new_dt = old_dt;
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
