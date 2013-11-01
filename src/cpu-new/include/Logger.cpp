#include "Logger.hpp"

Logger::Logger(int print_log, std::ofstream &out_file)
{

    gstream = (print_log ? &out_file : &std::cout);
}

Logger::~Logger()
{
}

void Logger::print_info(int n, double e2, double eta, float itime)
{
    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(5)  << std::left  << "N:";
    *gstream << std::setw(4)  << std::right << n;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(5)  << std::left  << "E2:";
    *gstream << std::setw(4)  << std::right << e2;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(5)  << std::left  << "Eta:";
    *gstream << std::setw(4)  << std::right << eta;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(5)  << std::left  << "T:";
    *gstream << std::setw(4)  << std::right << itime;
    *gstream << std::endl;

}

void Logger::print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f, double *dt)
{
    for (int i = 0; i < n; i++) {

            *gstream << std::fixed;
            gstream->precision(4);
            *gstream << std::setw(6) << ITIME;
            gstream->precision(2);
            *gstream << std::setw(6)  << i;

            // Scientific notation
            *gstream << std::scientific;
            gstream->precision(6);
            *gstream << std::setw(15) << std::right << r[i].w;

            *gstream << std::setw(15) << std::right << r[i].x;
            *gstream << std::setw(15) << std::right << r[i].y;
            *gstream << std::setw(15) << std::right << r[i].z;

            *gstream << std::setw(15) << std::right << v[i].x;
            *gstream << std::setw(15) << std::right << v[i].y;
            *gstream << std::setw(15) << std::right << v[i].z;

            *gstream << std::setw(15) << std::right << f[i].a[0];
            *gstream << std::setw(15) << std::right << f[i].a[1];
            *gstream << std::setw(15) << std::right << f[i].a[2];

            *gstream << std::setw(15) << std::right << f[i].a1[0];
            *gstream << std::setw(15) << std::right << f[i].a1[1];
            *gstream << std::setw(15) << std::right << f[i].a1[2];

            *gstream << std::setw(15) << std::right << dt[i];
            *gstream << std::endl;
    }
}

void Logger::print_energy_log(double ITIME, int iterations, long long interactions, int nsteps, Gtime gtime, Energy &energy, double new_energy)
{
    energy.end = new_energy;
    double rel_error     = std::abs((energy.end-energy.tmp)/energy.tmp);
    double cum_error     = std::abs((energy.end-energy.ini)/energy.ini);
    energy.tmp = energy.end;

    float time = omp_get_wtime() - gtime.integration_ini;
    if (ITIME != 0.0){
        gtime.gflops =  60.10e-9 * (interactions / gtime.update_end);
    }
    if(ITIME == 0.0)
    {
            *gstream << std::setw(2)  << std::left  << "#";
            *gstream << std::setw(13) << std::right << "IteTime";
            *gstream << std::setw(8)  << std::right << "Iter";
            *gstream << std::setw(10) << std::right << "Nsteps";
            *gstream << std::setw(16) << std::right << "ElapsedTime";

            *gstream << std::setw(15) << std::right << "Energy";
            *gstream << std::setw(15) << std::right << "RelE";
            *gstream << std::setw(15) << std::right << "CumE";

            *gstream << std::setw(12)  << std::right << "GFLOPS";
            *gstream << std::endl;
    }

    *gstream << std::fixed;
    gstream->precision(7);
    *gstream << std::setw(2)  << std::right << "00";
    *gstream << std::setw(13) << std::right << ITIME;
    *gstream << std::setw(8)  << std::right << iterations;
    *gstream << std::setw(10) << std::right << nsteps;
    *gstream << std::setw(16) << std::right << time;

    *gstream << std::scientific;
    gstream->precision(6);
    *gstream << std::setw(15) << std::right << energy.end;
    *gstream << std::setw(15) << std::right << rel_error;
    *gstream << std::setw(15) << std::right << cum_error;

    *gstream << std::fixed;
    gstream->precision(3);
    *gstream << std::setw(12) << std::right << gtime.gflops;
    *gstream << std::endl;
}
