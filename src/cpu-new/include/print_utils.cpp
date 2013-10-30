#include "print_utils.hpp"


void print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f, double *dt)
{
    for (int i = 0; i < n; i++) {

            std::cout << std::fixed;
            std::cout.precision(4);
            std::cout << std::setw(6) << ITIME;
            std::cout.precision(2);
            std::cout << std::setw(6)  << i;

            // Scientific notation
            std::cout << std::scientific;
            std::cout.precision(6);
            std::cout << std::setw(15) << std::right << r[i].w;
            std::cout << std::setw(15) << std::right << r[i].x;
            std::cout << std::setw(15) << std::right << r[i].y;
            std::cout << std::setw(15) << std::right << r[i].z;
            std::cout << std::setw(15) << std::right << v[i].x;
            std::cout << std::setw(15) << std::right << v[i].y;
            std::cout << std::setw(15) << std::right << v[i].z;
            std::cout << std::setw(15) << std::right << f[i].a[0];
            std::cout << std::setw(15) << std::right << f[i].a[1];
            std::cout << std::setw(15) << std::right << f[i].a[2];
            std::cout << std::setw(15) << std::right << f[i].a1[0];
            std::cout << std::setw(15) << std::right << f[i].a1[1];
            std::cout << std::setw(15) << std::right << f[i].a1[2];
            std::cout << std::setw(15) << std::right << dt[i];
            std::cout << std::endl;
    }
}


void print_energy_log(double ITIME, int iterations, long long interactions, int nsteps, Gtime gtime, Energy &energy, double new_energy)
{
    energy.end = new_energy;
    double rel_error     = std::abs((energy.end-energy.tmp)/energy.tmp);
    double rel_cum_error = std::abs((energy.end-energy.tmp)/energy.ini);
    double cum_error     = std::abs((energy.end-energy.ini)/energy.ini);
    energy.tmp = energy.end;

    float time = omp_get_wtime() - gtime.integration_ini;
    if (ITIME != 0.0){
        gtime.gflops =  60.10e-9 * (interactions / gtime.update_end);
    }
    if(ITIME == 0.0)
    {
            std::cout << std::setw(2)  << std::left  << "#";
            std::cout << std::setw(10) << std::right << "IteTime";
            std::cout << std::setw(8)  << std::right << "Iter";
            std::cout << std::setw(10) << std::right << "Nsteps";
            std::cout << std::setw(16) << std::right << "PredictionT";
            std::cout << std::setw(16) << std::right << "UpdateT";
            std::cout << std::setw(16) << std::right << "CorrectionT";
            std::cout << std::setw(16) << std::right << "ElapsedTime";
            std::cout << std::setw(12) << std::right << "Energy";
            std::cout << std::setw(15) << std::right << "RelE";
            std::cout << std::setw(15) << std::right << "CumRelE";
            std::cout << std::setw(15) << std::right << "CumE";
            std::cout << std::setw(12)  << std::right << "GFLOPS";
            std::cout << std::endl;
    }

        //printf("00 %7.3f %8d %10d %12.4f %12.4f %12.4f %12.4f %6.4f %5.3e %5.3e %5.3e %8.3f\n",

    std::cout << std::fixed;
    std::cout.precision(7);
    std::cout << std::setw(2)  << std::right << "00";
    std::cout << std::setw(10) << std::right << ITIME;
    std::cout << std::setw(8)  << std::right << iterations;
    std::cout << std::setw(10) << std::right << nsteps;
    std::cout << std::setw(16) << std::right << gtime.prediction_end;
    std::cout << std::setw(16) << std::right << gtime.update_end;
    std::cout << std::setw(16) << std::right << gtime.correction_end;
    std::cout << std::setw(16) << std::right << time;
    std::cout << std::setw(12) << std::right << energy.end;
    std::cout << std::scientific;
    std::cout.precision(6);
    std::cout << std::setw(15) << std::right << rel_error;
    std::cout << std::setw(15) << std::right << rel_cum_error;
    std::cout << std::setw(15) << std::right << cum_error;
    std::cout << std::fixed;
    std::cout.precision(3);
    std::cout << std::setw(12) << std::right << gtime.gflops;
    std::cout << std::endl;
}
