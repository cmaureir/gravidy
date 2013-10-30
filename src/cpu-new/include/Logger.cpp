#include "Logger.hpp"

Logger::Logger(int print_log, std::ofstream &out_file)
{

    gstream = (print_log ? &out_file : &std::cout);
}

Logger::~Logger()
{
    //*gstream.close();
}

void Logger::print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f, double *dt)
{
    for (int i = 0; i < n; i++) {

            //gstream << std::fixed;
            //gstream.precision(4);
            //gstream << std::setw(6) << ITIME;
            //gstream.precision(2);
            //gstream << std::setw(6)  << i;

            //// Scientific notation
            //gstream << std::scientific;
            //gstream.precision(6);
            //gstream << std::setw(15) << std::right << r[i].w;
            //gstream << std::setw(15) << std::right << r[i].x;
            //gstream << std::setw(15) << std::right << r[i].y;
            //gstream << std::setw(15) << std::right << r[i].z;
            //gstream << std::setw(15) << std::right << v[i].x;
            //gstream << std::setw(15) << std::right << v[i].y;
            //gstream << std::setw(15) << std::right << v[i].z;
            //gstream << std::setw(15) << std::right << f[i].a[0];
            //gstream << std::setw(15) << std::right << f[i].a[1];
            //gstream << std::setw(15) << std::right << f[i].a[2];
            //gstream << std::setw(15) << std::right << f[i].a1[0];
            //gstream << std::setw(15) << std::right << f[i].a1[1];
            //gstream << std::setw(15) << std::right << f[i].a1[2];
            //gstream << std::setw(15) << std::right << dt[i];
            //gstream << std::endl;
    }
}

void Logger::print_energy_log(double ITIME, int iterations, long long interactions, int nsteps, Gtime gtime, Energy &energy, double new_energy)
{
    *gstream << "hola\n";
//    energy.end = new_energy;
//    double rel_error     = std::abs((energy.end-energy.tmp)/energy.tmp);
//    double rel_cum_error = std::abs((energy.end-energy.tmp)/energy.ini);
//    double cum_error     = std::abs((energy.end-energy.ini)/energy.ini);
//    energy.tmp = energy.end;
//
//    float time = omp_get_wtime() - gtime.integration_ini;
//    if (ITIME != 0.0){
//        gtime.gflops =  60.10e-9 * (interactions / gtime.update_end);
//    }
//    if(ITIME == 0.0)
//    {
//            std::gstream << std::setw(2)  << std::left  << "#";
//            std::gstream << std::setw(10) << std::right << "IteTime";
//            std::gstream << std::setw(8)  << std::right << "Iter";
//            std::gstream << std::setw(10) << std::right << "Nsteps";
//            std::gstream << std::setw(16) << std::right << "PredictionT";
//            std::gstream << std::setw(16) << std::right << "UpdateT";
//            std::gstream << std::setw(16) << std::right << "CorrectionT";
//            std::gstream << std::setw(16) << std::right << "ElapsedTime";
//            std::gstream << std::setw(12) << std::right << "Energy";
//            std::gstream << std::setw(15) << std::right << "RelE";
//            std::gstream << std::setw(15) << std::right << "CumRelE";
//            std::gstream << std::setw(15) << std::right << "CumE";
//            std::gstream << std::setw(12)  << std::right << "GFLOPS";
//            std::gstream << std::endl;
//    }
//
//    std::gstream << std::fixed;
//    std::gstream.precision(7);
//    std::gstream << std::setw(2)  << std::right << "00";
//    std::gstream << std::setw(10) << std::right << ITIME;
//    std::gstream << std::setw(8)  << std::right << iterations;
//    std::gstream << std::setw(10) << std::right << nsteps;
//    std::gstream << std::setw(16) << std::right << gtime.prediction_end;
//    std::gstream << std::setw(16) << std::right << gtime.update_end;
//    std::gstream << std::setw(16) << std::right << gtime.correction_end;
//    std::gstream << std::setw(16) << std::right << time;
//    std::gstream << std::setw(12) << std::right << energy.end;
//    std::gstream << std::scientific;
//    std::gstream.precision(6);
//    std::gstream << std::setw(15) << std::right << rel_error;
//    std::gstream << std::setw(15) << std::right << rel_cum_error;
//    std::gstream << std::setw(15) << std::right << cum_error;
//    std::gstream << std::fixed;
//    std::gstream.precision(3);
//    std::gstream << std::setw(12) << std::right << gtime.gflops;
//    std::gstream << std::endl;
}
