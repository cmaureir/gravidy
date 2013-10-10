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
    for (int i = 0; i < n; i++)
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
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double a2 = get_magnitude(h_f[i].a[0],  h_f[i].a[1],  h_f[i].a[2]);
        double j2 = get_magnitude(h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2]);
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


/*
 * @fn save_old()
 *
 * @brief
 *  Save the previous particle position, velocity,
 *    acceleration and jrk.
 */
void save_old(int total)
{
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
    gtime.prediction_ini = omp_get_wtime();
    //#pragma omp parallel for
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double dt  = ITIME - h_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        h_p[i].r[0] = (dt3/6 * h_f[i].a1[0]) + (dt2/2 * h_f[i].a[0]) + (dt * h_v[i].x) + h_r[i].x;
        h_p[i].r[1] = (dt3/6 * h_f[i].a1[1]) + (dt2/2 * h_f[i].a[1]) + (dt * h_v[i].y) + h_r[i].y;
        h_p[i].r[2] = (dt3/6 * h_f[i].a1[2]) + (dt2/2 * h_f[i].a[2]) + (dt * h_v[i].z) + h_r[i].z;

        h_p[i].v[0] = (dt2/2 * h_f[i].a1[0]) + (dt * h_f[i].a[0]) + h_v[i].x;
        h_p[i].v[1] = (dt2/2 * h_f[i].a1[1]) + (dt * h_f[i].a[1]) + h_v[i].y;
        h_p[i].v[2] = (dt2/2 * h_f[i].a1[2]) + (dt * h_f[i].a[2]) + h_v[i].z;

        h_p[i].m = h_m[i];

    }
    gtime.prediction_end += omp_get_wtime() - gtime.prediction_ini;
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
        h_a2[i].x = (-6 * (h_old_a[i].x - h_f[i].a[0] ) - dt1 * (4 * h_old_a1[i].x + 2 * h_f[i].a1[0]) ) / dt2;
        h_a2[i].y = (-6 * (h_old_a[i].y - h_f[i].a[1] ) - dt1 * (4 * h_old_a1[i].y + 2 * h_f[i].a1[1]) ) / dt2;
        h_a2[i].z = (-6 * (h_old_a[i].z - h_f[i].a[2] ) - dt1 * (4 * h_old_a1[i].z + 2 * h_f[i].a1[2]) ) / dt2;

        // Acceleration 3rd derivate
        h_a3[i].x = (12 * (h_old_a[i].x - h_f[i].a[0] ) + 6 * dt1 * (h_old_a1[i].x + h_f[i].a1[0]) ) / dt3;
        h_a3[i].y = (12 * (h_old_a[i].y - h_f[i].a[1] ) + 6 * dt1 * (h_old_a1[i].y + h_f[i].a1[1]) ) / dt3;
        h_a3[i].z = (12 * (h_old_a[i].z - h_f[i].a[2] ) + 6 * dt1 * (h_old_a1[i].z + h_f[i].a1[2]) ) / dt3;

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
        normal_dt = normalize_dt(normal_dt, h_dt[i], h_t[i], i);
        h_dt[i] = normal_dt;

    }

    // Copy the new r and v to from the device to the host
    CUDA_SAFE_CALL(cudaMemcpy(d_r, h_r, d4_size, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_v, h_v, d4_size, cudaMemcpyHostToDevice));

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

/*
 * @fn normalize_dt()
 *
 * @brief
 *  Normalize the previous value of the time step of a particle i,
 *  comparing it with its previous value and fitting the new value
 *  in the blocks before of after the current block.
 *
 */
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
