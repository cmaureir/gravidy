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
            printf("%.10f %6d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
                    ITIME,
                    i,
                    h_m[i],
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
            fprintf(out,"%.10f\t%6d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
                    ITIME,
                    i,
                    h_m[i],
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
    printf("% 5d % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f\n",
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
 *        * Total iterations
 *        * Amount of particles which were updated
 *        * Prediction time
 *        * Update time
 *        * Correction time
 *        * Total clock time
 *        * Energy
 *        * Relative Energy Error      (E_{t} - E_{t-1})/E_{t-1}
 *        * Relative Cumulative Error  (E_{t} - E_{t-1})/E_{0}
 *        * Cumulative Energy Error    (E_{t} - E_{0})/E_{0}
 *        * GFLOPS of the k_update kernel only.
 *
 */
void get_energy_log(double ITIME, int iterations, long long interactions, int nsteps, FILE *out, double energy)
{
    energy_end = energy;
    double rel_error     = abs((energy_end-energy_tmp)/energy_tmp);
    double rel_cum_error = abs((energy_end-energy_tmp)/energy_ini);
    double cum_error     = abs((energy_end-energy_ini)/energy_ini);
    energy_tmp = energy_end;
    float time = omp_get_wtime() - gtime.integration_ini;
    if (ITIME != 0.0){
        gflops =  60.10e-9 * (interactions / gtime.update_end);
    }
    if(ITIME == 0.0)
    {
        if (print_log)
        {
            fprintf(out,"00 %7s %8s %10s %12s %12s %12s %12s %7s %9s %9s %9s %8s\n",
                "Time",
                "Iterations",
                "Nsteps",
                "PredictionT",
                "UpdateT",
                "CorrectionT",
                "TTime",
                "Energy",
                "RelErr",
                "CumRelErr",
                "CumErr",
                "GFLOPS");
            fflush(out);
        }
        else
        {

            printf("00 %7s %8s %10s %12s %12s %12s %12s %7s %9s %9s %9s %8s\n",
                      "IteTime",
                      "Iter",
                      "Nsteps",
                      "PredictionT",
                      "UpdateT",
                      "CorrectionT",
                      "ElapsedTime",
                      "Energy",
                      "RelE",
                      "CumRelE",
                      "CumE",
                      "GFLOPS");
        }
    }

    if (print_log)
    {
        fprintf(out,"00 %7.3f %8d %10d %12.4f %12.4f %12.4f %12.4f %6.4f %5.3e %5.3e %5.3e %8.3f\n",
            ITIME,
            iterations,
            nsteps,
            gtime.prediction_end,
            gtime.update_end,
            gtime.correction_end,
            time,
            energy_end,
            rel_error,
            rel_cum_error,
            cum_error,
            gflops);
        print_all(n,ITIME,out);
        fflush(out);
    }
    else
    {
        printf("00 %7.3f %8d %10d %12.4f %12.4f %12.4f %12.4f %6.4f %5.3e %5.3e %5.3e %8.3f\n",
            ITIME,
            iterations,
            nsteps,
            gtime.prediction_end,
            gtime.update_end,
            gtime.correction_end,
            time,
            energy_end,
            rel_error,
            rel_cum_error,
            cum_error,
            gflops);
    }
}

string get_time(){

    time_t timeObj;
    char buffer[100];

    time(&timeObj);
    tm *pTime = localtime(&timeObj);
    sprintf(buffer, "%02d-%02d-%04d_%02d:%02d:%02d", pTime->tm_mday,
                                                     pTime->tm_mon+1,
                                                     1900+pTime->tm_year,
                                                     pTime->tm_hour,
                                                     pTime->tm_min,
                                                     pTime->tm_sec);
    return buffer;
}
