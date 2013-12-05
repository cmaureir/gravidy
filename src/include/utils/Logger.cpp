#include "Logger.hpp"

Logger::Logger(NbodySystem *ns)
{
    this->ns = ns;
    gstream = (this->ns->ops.print_screen ? &std::cout : &this->ns->out_file);
}

Logger::~Logger()
{
    // None
}

void Logger::print_info()
{
    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(8)  << std::left  << "N:";
    *gstream << std::setw(8)  << std::right << ns->n;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(8)  << std::left  << "E2:";
    *gstream << std::setw(8)  << std::right << ns->e2;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(8)  << std::left  << "Eta:";
    *gstream << std::setw(8)  << std::right << ns->eta;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(8)  << std::left  << "T:";
    *gstream << std::setw(8)  << std::right << ns->integration_time;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(8)  << std::left  << "T_rh:";
    *gstream << std::setw(8)  << std::right << ns->hmr_time;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(8)  << std::left  << "T_cr:";
    *gstream << std::setw(8)  << std::right << ns->cr_time;
    *gstream << std::endl;


}

void Logger::print_lagrange_radii(double ITIME, std::vector<double> lagrange_radii)
{
    *gstream << "01 ";
    *gstream << std::fixed;
    *gstream << ITIME << " ";
    gstream->precision(4);
    for (int i = 0; i < 1/RADIUS_RATIO; i++) {
        *gstream << std::setw(6) << std::right << lagrange_radii[i] << " ";
    }
    *gstream << std::endl;
}

void Logger::print_all(double ITIME)
{
    for (int i = 0; i < ns->n; i++) {

            *gstream << std::fixed;
            gstream->precision(4);
            *gstream << std::setw(6) << ITIME;
            gstream->precision(2);
            *gstream << std::setw(6)  << i;

            // Scientific notation
            *gstream << std::scientific;
            gstream->precision(6);
            *gstream << std::setw(15) << std::right << ns->h_r[i].w;

            *gstream << std::setw(15) << std::right << ns->h_r[i].x;
            *gstream << std::setw(15) << std::right << ns->h_r[i].y;
            *gstream << std::setw(15) << std::right << ns->h_r[i].z;

            *gstream << std::setw(15) << std::right << ns->h_v[i].x;
            *gstream << std::setw(15) << std::right << ns->h_v[i].y;
            *gstream << std::setw(15) << std::right << ns->h_v[i].z;

            *gstream << std::setw(15) << std::right << ns->h_f[i].a[0];
            *gstream << std::setw(15) << std::right << ns->h_f[i].a[1];
            *gstream << std::setw(15) << std::right << ns->h_f[i].a[2];

            *gstream << std::setw(15) << std::right << ns->h_f[i].a1[0];
            *gstream << std::setw(15) << std::right << ns->h_f[i].a1[1];
            *gstream << std::setw(15) << std::right << ns->h_f[i].a1[2];

            *gstream << std::setw(15) << std::right << ns->h_dt[i];
            *gstream << std::endl;
    }
}

void Logger::print_energy_log(double ITIME, int iterations, long long interactions, int nsteps, double new_energy)
{
    ns->en.end = new_energy;
    double rel_error     = std::abs((ns->en.end - ns->en.tmp)/ns->en.tmp);
    double cum_error     = std::abs((ns->en.end - ns->en.ini)/ns->en.ini);
    ns->en.tmp = ns->en.end;

    float time = omp_get_wtime() - ns->gtime.integration_ini;
    if (ITIME != 0.0){
        ns->gtime.gflops =  60.10e-9 * (interactions / ns->gtime.update_end);
    }
    if(ITIME == 0.0)
    {
            *gstream << std::setw(2)  << std::left  << "#";
            *gstream << std::setw(10) << std::right << "IteTime";
            *gstream << std::setw(15) << std::right << "Iter";
            *gstream << std::setw(15) << std::right << "Nsteps";

            *gstream << std::setw(20) << std::right << "Energy";
            *gstream << std::setw(15) << std::right << "RelE";
            *gstream << std::setw(15) << std::right << "CumE";
            *gstream << std::setw(16) << std::right << "ElapsedTime";

            *gstream << std::setw(12)  << std::right << "GFLOPS";
            *gstream << std::endl;
    }

    *gstream << std::fixed;
    *gstream << std::setw(2)  << std::right << "00";
    gstream->precision(3);
    *gstream << std::setw(10) << std::right << ITIME;
    *gstream << std::setw(15)  << std::right << iterations;
    *gstream << std::setw(15) << std::right << nsteps;

    *gstream << std::scientific;
    gstream->precision(10);
    *gstream << std::setw(20) << std::right << ns->en.end;
    gstream->precision(6);
    *gstream << std::setw(15) << std::right << rel_error;
    *gstream << std::setw(15) << std::right << cum_error;
    *gstream << std::setw(16) << std::right << time;

    *gstream << std::fixed;
    gstream->precision(3);
    *gstream << std::setw(12) << std::right << ns->gtime.gflops;
    *gstream << std::endl;
}
