#include "dynamics_cpu.hpp"

/*
 * @fn next_itime()
 *
 * @brief
 *  Determine the next iteration time
 *  considering the time of all the particles
 */
void
next_itime(float *ATIME)
{
    *ATIME = 1.0e10;
    for (int i = 0; i < n; i++)
    {
        *ATIME = std::min(*ATIME,h_t[i] + h_dt[i]);
    }
}

/*
 * @fn find_particles_to_move()
 *
 * @brief
 *  Find particles with the same time that the
 *  integration time.
 *
 * @todo
 *  try if is better to kept the index in the same position
 */
int find_particles_to_move(float ITIME)
{
    int i, j;
    j = 0;
    for (i = 0; i < n; i++)
    {
        h_move[i] = -1;
        if (h_t[i] + h_dt[i] == ITIME)
        {
            h_move[j] = i;
            j++;
//            std::cout << i << " ";
        }
    }
//    std::cout << std::endl;
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
 *
 */
void
init_dt(float *ATIME)
{
    int i;
    float tmp_dt;

    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    for (i = 0; i < n; i++) {
        tmp_dt = ETA_S *
                 magnitude(h_a[i].x, h_a[i].y, h_a[i].z) /
                 magnitude(h_j[i].x, h_j[i].y, h_j[i].z);

        /* Adjusting to block timesteps */
        if (tmp_dt < D_TIME_MIN)
        {
            tmp_dt = D_TIME_MIN;
        }
        else if (tmp_dt > D_TIME_MAX)
        {
            tmp_dt = D_TIME_MAX;
        }
        else
        {
          tmp_dt = pow(2,(int)((log(tmp_dt)/log(2.0))-1));
        }
        *ATIME = std::min(*ATIME, tmp_dt);
        h_dt[i] = tmp_dt;
        h_t[i] = 0.0;
    }

}


/*
 * @fn init_acc_jerk()
 *
 * @brief
 *  Initial calculation of the all the acceleration and jerk.
 */
void
init_acc_jerk()
{
    double rx, ry, rz;
    double vx, vy, vz;
    double f, f3, f5;
    int i, j;


    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            rx = h_r[j].x - h_r[i].x;
            ry = h_r[j].y - h_r[i].y;
            rz = h_r[j].z - h_r[i].z;

            vx = h_v[j].x - h_v[i].x;
            vy = h_v[j].y - h_v[i].y;
            vz = h_v[j].z - h_v[i].z;

            f = rx*rx + ry*ry + rz*rz + E*E;

            f3 = f * f * f;
            f5 = f3 * f * f;
            f3 = sqrt(f3);
            f5 = sqrt(f5);

            h_a[i].x += h_m[j] * rx / f3;
            h_a[i].y += h_m[j] * ry / f3;
            h_a[i].z += h_m[j] * rz / f3;

            h_j[i].x += h_m[j] * (vx/f3 + (3 * vx * rx * rx)/f5);
            h_j[i].y += h_m[j] * (vy/f3 + (3 * vy * ry * ry)/f5);
            h_j[i].z += h_m[j] * (vz/f3 + (3 * vz * rz * rz)/f5);
        }
    }
}


/*
 * @fn update_acc_jerk()
 *
 * @brief
 *  Iteration calculation of the acceleration and the jerk.
 *
 * @todo
 *  Only perform the calculation of the particles with the t + dt = ITIME
 *
 */
void
update_acc_jerk(int total)
{
    double rx, ry, rz;
    double vx, vy, vz;
    double f, f3, f5;
    int i, j, k;


    for (k = 0; k < total; k++)
    {
        i = h_move[k];
        // Cleaning acceleration and jerk
        h_a[i].x = h_a[i].y = h_a[i].z = 0.0f;
        h_j[i].x = h_j[i].y = h_j[i].z = 0.0f;

        for (j = 0; j < n; j++)
        {
            if(i != j)
            {
                rx = h_p_r[j].x - h_p_r[i].x;
                ry = h_p_r[j].y - h_p_r[i].y;
                rz = h_p_r[j].z - h_p_r[i].z;

                vx = h_p_v[j].x - h_p_v[i].x;
                vy = h_p_v[j].y - h_p_v[i].y;
                vz = h_p_v[j].z - h_p_v[i].z;

                if(i == 0 || j == 0 || i == 1 || j == 1)
                {
                    #define NEW_E (1e-8)
                    f = rx*rx + ry*ry + rz*rz + NEW_E*NEW_E;
                }
                else
                {
                    f = rx*rx + ry*ry + rz*rz + E*E;
                }

                f3 = f * f * f;
                f5 = f3 * f * f;
                f3 = sqrt(f3);
                f5 = sqrt(f5);

                h_a[i].x += h_m[j] * rx / f3;
                h_a[i].y += h_m[j] * ry / f3;
                h_a[i].z += h_m[j] * rz / f3;

                h_j[i].x += h_m[j] * (vx/f3 + (3 * vx * rx * rx)/f5);
                h_j[i].y += h_m[j] * (vy/f3 + (3 * vy * ry * ry)/f5);
                h_j[i].z += h_m[j] * (vz/f3 + (3 * vz * rz * rz)/f5);
            }
        }

    }
}

/*
 * @fn energy()
 *
 * @brief
 *  Kinetic energy calculation.
 *
 *  The Potential energy is calculated in the init_acc_jerk()
 *    and in the update_acc_jerk() functions.
 *
 */
float
energy()
{
    float ekin_tmp = 0;
    float epot_tmp = 0;
    double r2, v2;
    double rx, ry, rz;
    double vx, vy, vz;
    int i, j;

    epot = 0;
    ekin = 0;

    for (i = 0; i < n; i++)
    {
        epot_tmp = 0;
        for (j = i+1; j < n; j++)
        {
            rx = h_p_r[j].x - h_p_r[i].x;
            ry = h_p_r[j].y - h_p_r[i].y;
            rz = h_p_r[j].z - h_p_r[i].z;

            r2 = rx*rx + ry*ry + rz*rz;// + E*E;

            epot_tmp -= (h_m[i] * h_m[j]) / sqrt(r2);
        }

        vx = h_p_v[i].x * h_p_v[i].x;
        vy = h_p_v[i].y * h_p_v[i].y;
        vz = h_p_v[i].z * h_p_v[i].z;

        v2 = vx + vy + vz;

        ekin_tmp = 0.5 * h_m[i] * v2;
        ekin += ekin_tmp;
        epot += epot_tmp;
    }

    return epot + ekin;
}

/*
 * @fn initial_energy()
 *
 * @brief
 *  Kinetic energy calculation.
 *
 *  The Potential energy is calculated in the init_acc_jerk()
 *    and in the update_acc_jerk() functions.
 *
 */
float
initial_energy()
{
    float ekin_tmp = 0;
    float epot_tmp = 0;
    double r2, v2;
    double rx, ry, rz;
    double vx, vy, vz;
    int i, j;

    epot = 0;
    ekin = 0;

    for (i = 0; i < n; i++)
    {
        epot_tmp = 0;
        for (j = i+1; j < n; j++)
        {
            rx = h_r[j].x - h_r[i].x;
            ry = h_r[j].y - h_r[i].y;
            rz = h_r[j].z - h_r[i].z;

            r2 = rx*rx + ry*ry + rz*rz;// + E*E;

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

void get_energy_log(int OUT, float ITIME)
{
    if(iterations % OUT == 0)
    {
        fprintf(stderr, "%.10f %.10e\n", ITIME, energy_total/OUT);
        print_positions(100);
        energy_total = 0.0f;
        energy_ini = energy_end;
    }
}

/*
 * @fn save_old()
 *
 * @brief
 *  Save the previous particle position, velocity,
 *    acceleration and jerk.
 *
 */
void
save_old()
{
    //#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
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
predicted_pos_vel(float ITIME)
{
    float tmp_dt, tmp_dt2, tmp_dt3;

    for (int i = 0; i < n; i++)
    {
        tmp_dt = ITIME - h_t[i];
        tmp_dt2 = (tmp_dt  * tmp_dt)/2;
        tmp_dt3 = (tmp_dt2 * tmp_dt)/6;

        h_p_r[i].x = tmp_dt3 * h_j[i].x + tmp_dt2 * h_a[i].x + tmp_dt * h_v[i].x + h_r[i].x;
        h_p_r[i].y = tmp_dt3 * h_j[i].y + tmp_dt2 * h_a[i].y + tmp_dt * h_v[i].y + h_r[i].y;
        h_p_r[i].z = tmp_dt3 * h_j[i].z + tmp_dt2 * h_a[i].z + tmp_dt * h_v[i].z + h_r[i].z;

        h_p_v[i].x = tmp_dt2 * h_j[i].x + tmp_dt * h_a[i].x + h_v[i].x;
        h_p_v[i].y = tmp_dt2 * h_j[i].y + tmp_dt * h_a[i].y + h_v[i].y;
        h_p_v[i].z = tmp_dt2 * h_j[i].z + tmp_dt * h_a[i].z + h_v[i].z;
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
void
predicted_pos_vel_kepler(float ITIME)
{
    // TO DO
}

/*
 * @fn correction_pos_vel()
 *
 * @brief
 *  Perform the calculation of the acceleration derivates
 *    to apply the correction terms to the position and velocity.
 *
 */
void
correction_pos_vel(float ITIME, int total)
{

    // Acceleration derivates (2nd and 3rd)
    double ax0_2, ay0_2, az0_2;
    double ax0_3, ay0_3, az0_3;
    double ax1_2, ay1_2, az1_2;

    float dt1, dt2, dt3, dt4, dt5;
    double abs_a1, abs_j1, abs_j12, abs_a1_3, abs_a1_2, abs_a1_22;

    float tmp_dt;
    int i, k;


    //#pragma omp parallel for
    for (k = 0; k < total; k++)
    {
        i = h_move[k];

        dt1 = ITIME - h_t[i];
        h_t[i] = ITIME;
        dt2 = dt1 * dt1;
        dt3 = dt2 * dt1;
        dt4 = dt2 * dt2;
        dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        ax0_2 = (-6 * (h_old_a[i].x - h_a[i].x ) - dt1 * (4 * h_old_j[i].x + 2 * h_j[i].x) ) / dt2;
        ay0_2 = (-6 * (h_old_a[i].y - h_a[i].y ) - dt1 * (4 * h_old_j[i].y + 2 * h_j[i].y) ) / dt2;
        az0_2 = (-6 * (h_old_a[i].z - h_a[i].z ) - dt1 * (4 * h_old_j[i].z + 2 * h_j[i].z) ) / dt2;

        // Acceleration 3rd derivate
        ax0_3 = (12 * (h_old_a[i].x - h_a[i].x ) + 6 * dt1 * (h_old_j[i].x + h_j[i].x) ) / dt3;
        ay0_3 = (12 * (h_old_a[i].y - h_a[i].y ) + 6 * dt1 * (h_old_j[i].y + h_j[i].y) ) / dt3;
        az0_3 = (12 * (h_old_a[i].z - h_a[i].z ) + 6 * dt1 * (h_old_j[i].z + h_j[i].z) ) / dt3;

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
        ax1_2 = ax0_2 + dt1 * ax0_3;
        ay1_2 = ay0_2 + dt1 * ay0_3;
        az1_2 = az0_2 + dt1 * az0_3;

        // |a_{1,i}|
        abs_a1   = sqrt((h_a[i].x * h_a[i].x) + (h_a[i].y * h_a[i].y) + (h_a[i].z * h_a[i].z));
        // |j_{1,i}|
        abs_j1   = sqrt((h_j[i].x * h_j[i].x) + (h_j[i].y * h_j[i].y) + (h_j[i].z * h_j[i].z));
        // |j_{1,i}|^{2}
        abs_j12  = abs_j1 * abs_j1;
        // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
        abs_a1_3 = sqrt((ax0_3 * ax0_3) + (ay0_3 * ay0_3) + (az0_3 * az0_3));
        // |a_{1,i}^{(2)}|
        abs_a1_2 = sqrt((ax1_2 * ax1_2) + (ay1_2 * ay1_2) + (az1_2 * az1_2));
        // |a_{1,i}^{(2)}|^{2}
        abs_a1_22  = abs_a1_2 * abs_a1_2;

        tmp_dt = sqrt( ETA_N * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));

        /* Adjusting to block timesteps */
        if (tmp_dt < D_TIME_MIN)
        {
            tmp_dt = D_TIME_MIN;
        }
        else if (tmp_dt > D_TIME_MAX)
        {
            tmp_dt = D_TIME_MAX;
        }
        else
        {
            tmp_dt = pow(2,(int)((log(tmp_dt)/log(2.0))-1));
        }
        h_dt[i] = tmp_dt;
    }
}


double magnitude(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}


