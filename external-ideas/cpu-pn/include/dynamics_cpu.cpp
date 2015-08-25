#include "dynamics_cpu.hpp"

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
    int i,j;
    #pragma omp parallel for private(j) schedule(dynamic, 24)
    for (i = INIT_PARTICLE; i < n; i++)
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
	gtime.update_ini = omp_get_wtime();
    int i, j, k;
    #pragma omp parallel for private(i,j)
    for (k = 0; k < total; k++)
    {
        i = h_move[k];
        h_f[i].a[0]  = 0.0;
        h_f[i].a[1]  = 0.0;
        h_f[i].a[2]  = 0.0;
        h_f[i].a1[0] = 0.0;
        h_f[i].a1[1] = 0.0;
        h_f[i].a1[2] = 0.0;

        //#pragma omp parallel for
        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i,j);
        }
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


    gtime.correction_end += omp_get_wtime() - gtime.correction_ini;
}
