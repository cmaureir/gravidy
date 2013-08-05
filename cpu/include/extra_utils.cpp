#include "extra_utils.hpp"

/*
 * @fn print_all()
 *
 * @brief Print all the information of each Particle.
 *
 */
void print_all(int limit, float ITIME, FILE *out)
{

    if (out == NULL)
    {
        for (int i = 0; i < limit; i++) {
            printf("%.10f %6d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
                    ITIME,
                    i,
                    h_r[i].x, h_r[i].y, h_r[i].z,
                    h_v[i].x, h_v[i].y, h_v[i].z,
                    h_f[i].a[0], h_f[i].a[1], h_f[i].a[2],
                    h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2],
                    h_dt[i]);
        }
    }
    else if (out != NULL)
    {
        for (int i = 0; i < limit; i++) {
            fprintf(out,"%.10f %6d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
                    ITIME,
                    i,
                    h_r[i].x, h_r[i].y, h_r[i].z,
                    h_v[i].x, h_v[i].y, h_v[i].z,
                    h_f[i].a[0], h_f[i].a[1], h_f[i].a[2],
                    h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2],
                    h_dt[i]);
        }
    }
}

/*
 * @fn print_forces()
 *
 * @brief Print acceleration and its first derivative for all the particles
 *
 */
void print_forces(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.10f %.10f %.10f %.10f %.10f %.10f\n",
                i,
                h_f[i].a[0], h_f[i].a[1], h_f[i].a[2],
                h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2] );
    }

}

/*
 * @fn print_times()
 *
 * @brief Print the iteration time, timestep, power of 2 exponent and time
 *       for each particle of the system.
 *
 */
void print_times(int limit, float itime)
{
    int exp = 0;
    for (int i = 0; i < limit; i++) {
        exp = (int)std::ceil(log(h_dt[i])/log(2));
        printf("%.10f %6d %.10f %2d %.10f\n",
               itime, i, h_dt[i], exp, h_t[i]);
    }
}

/*
 * @fn print_predicted()
 *
 * @brief Print predicted position and velocity of all the particles
 *
 */
void print_predicted(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.10f %.10f %.10f %.10f %.10f %.10f\n",
                i,
                h_p[i].r[0], h_p[i].r[1], h_p[i].r[2],
                h_p[i].v[0], h_p[i].v[1], h_p[i].v[2]);
    }
}

/*
 * @fn print_movement
 *
 * @brief Print the particles which need to be updated in this iteration time
 *
 */
void print_movement(int limit, int total, float ITIME)
{
    for (int i = 0; i < limit; i++)
    {
        int value = 0;
        for(int j = 0; j < total ; j++)
        {
            int k = h_move[j];
            if (k == i)
                value = 1;

        }
        printf("%d ", value);
    }
    printf("\n");
}

/*
 * @fn print_particle()
 *
 * @brief Print information of a i-particle
 *
 */
void print_particle(int i)
{
    printf("%5d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
            i,
            h_r[i].x,     h_r[i].y,     h_r[i].z,
            h_v[i].x,     h_v[i].y,     h_v[i].z,
            h_f[i].a[0],  h_f[i].a[1],  h_f[i].a[2],
            h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2],
            h_dt[i]);

}


/*
 * @fn get_magnitude()
 *
 * @brief Calculate the magnitude of a 3D vector.
 *
 */
double get_magnitude(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

/*
 * @fn get_energy_log()
 *
 * @brief Print a snapshoot with information of the system.
 *        * N-body time
 *        * CPU iterations
 *        * GPU iterations
 *        * Total iterations
 *        * Amount of particles which were updated
 *        * GPU clock time
 *        * CPU clock time
 *        * Total clock time
 *        * Energy
 *        * Relative Energy Error      (E_{t} - E_{t-1})/E_{t-1}
 *        * Relative Cumulative Error  (E_{t} - E_{t-1})/E_{0}
 *        * Cumulative Energy Error    (E_{t} - E_{0})/E_{0}
 *
 */
void get_energy_log(double ITIME, int iterations, int nsteps, FILE *out, double energy)
{
    energy_end = energy;
    double rel_error     = abs((energy_end-energy_tmp)/energy_tmp);
    double rel_cum_error = abs((energy_end-energy_tmp)/energy_ini);
    double cum_error     = abs((energy_end-energy_ini)/energy_ini);
    energy_tmp = energy_end;
    float time = omp_get_wtime() - ini_time;

    if(ITIME == 0.0)
    {
       fprintf(out, "00 %3s\t %10s\t %10s\t %10s\t %10s\t %8s\t %8s\t %8s\t %8s\t %8s\t %8s\n",
     // printf(     "00 %3s\t %10s\t %10s\t %10s\t %10s\t %8s\t %8s\t %8s\t %8s\t %8s\t %8s\n",
                "Time",
                "CpuIte",
                "GpuIte",
                "Ite",
                "Nsteps",
                "GpuT",
                "TTime",
                "Energy",
                "RelErr",
                "CumRelErr",
                "CumErr");
    }
   fprintf(out, "00 % 3f\t % 10d\t % 10d\t % 10d\t % 10d\t % 6.4f\t % 6.4f\t % .6e\t % .6e\t % .6e\t % .6e\n",
//   printf(     "00 % 3f\t % 10d\t % 10d\t % 10d\t % 10d\t % 6.4f\t % 6.4f\t % .6e\t % .6e\t % .6e\t % .6e\n",
            ITIME,
            iterations,
            0,
            iterations,
            nsteps,
            0.0,
            time,
            energy_end,
            rel_error,
            rel_cum_error,
            cum_error);

    if (print_log)
    {
        print_all(n,ITIME,out);
    }
    fflush(out);
}

