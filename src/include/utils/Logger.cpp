#include "Logger.hpp"

Logger::Logger(NbodySystem *ns)
{
    this->ns = ns;
    this->print_screen = this->ns->ops.print_screen;
    this->ofname = "";

    // Stripping the output filename if it's inside a directory
    std::vector<std::string> tmp;
    std::stringstream ss(ns->output_filename);
    std::string item;
    while (std::getline(ss, item, '/')) {
        tmp.push_back(item);
    }
    std::string tmp_fname = this->get_timestamp() +"_"+ tmp[tmp.size()-1];

    for(int i = 0; i < (int)tmp.size() - 1;i++)
    {
        this->ofname += tmp[i] + "/";
    }
    this->ofname += tmp_fname;

}

Logger::~Logger()
{
    // None
}

std::string Logger::get_timestamp()
{
    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::ostringstream s;

    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_hour << ":";
    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_min  << ":";
    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_sec  << "_";
    s << std::setw(2) << std::setfill('0') << ltm->tm_mday << "-";
    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_mon<< "-";
    s << std::setw(4) << std::setfill('0') << 1900 + ltm->tm_year;

    return s.str();
}

void Logger::print_info()
{

    if(print_screen)
    {
        gstream = &std::cout;
    }
    else
    {
        std::string ofname_info = ofname + ".info";

        out_file.open(ofname_info.c_str(), std::ios::out);
        if(!out_file)
        {
          std::cerr << "gravidy: cannot open file "
                    << ofname_info.c_str()
                    << ": No such file or directory"
                    << std::endl;
          exit(1);
        }
        gstream = &out_file;
    }

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35)  << std::left  << "Number of particles:";
    *gstream << std::setw(15)  << std::right << ns->n;
    *gstream << std::endl;

    *gstream << std::scientific;
    *gstream << std::setprecision(4);

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35)  << std::left  << "Softening:";
    *gstream << std::setw(15)  << std::right << sqrt(ns->e2);
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35) << std::left  << "Eta (for timesteps):";
    *gstream << std::setw(15)  << std::right << ns->eta;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35)  << std::left  << "Integration time:";
    *gstream << std::setw(15)  << std::right << ns->integration_time;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35) << std::left  << "Half-mass radius (r_h):";
    *gstream << std::setw(15)  << std::right << ns->r_h;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35) << std::left  << "Center of Density (cod):";
    *gstream << std::setw(15)  << std::right << ns->cod.x << " " << ns->cod.y << " " << ns->cod.z;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35) << std::left  << "Half-mass relaxation time (T_rh):";
    *gstream << std::setw(15)  << std::right << ns->hmr_time;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35) << std::left  << "Half-mass relaxation time (soft):";
    *gstream << std::setw(15)  << std::right << ns->hmr_time_soft;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(35) << std::left  << "Crossing time (T_cr):";
    *gstream << std::setw(15)  << std::right << ns->cr_time;
    *gstream << std::endl;

    *gstream << std::fixed;

    if(!print_screen)
    {
        out_file.close();
    }

}

void Logger::print_lagrange_radii(double ITIME, std::vector<double> lagrange_radii)
{
    if(print_screen)
    {
        gstream = &std::cout;
    }
    else
    {
        std::string ofname_radii = ofname + ".radii";
        out_file.open(ofname_radii.c_str(), std::ios::out | std::ios::app );
        if(!out_file)
        {
          std::cerr << "gravidy: cannot open file "
                    << ofname_radii.c_str()
                    << ": No such file or directory"
                    << std::endl;
          exit(1);
        }
        gstream = &out_file;
    }


    *gstream << "01 ";
    *gstream << std::fixed;
    *gstream << ITIME << " ";
    gstream->precision(4);
    for (int i = 0; i < 1/RADIUS_RATIO - 1; i++) {
        *gstream << std::setw(6) << std::right << lagrange_radii[i] << " ";
    }
    *gstream << std::endl;

    if(!print_screen)
    {
        out_file.close();
    }
}

void Logger::print_all(double ITIME)
{
    if(print_screen)
    {
        gstream = &std::cout;
    }
    else
    {
        std::ostringstream s;
        s << std::setw(4) << std::setfill('0') << ITIME;
        std::string ofname_all = ofname + ".all.t" + s.str();
        out_file.open(ofname_all.c_str(), std::ios::out);
        if(!out_file)
        {
          std::cerr << "gravidy: cannot open file "
                    << ofname_all.c_str()
                    << ": No such file or directory"
                    << std::endl;
          exit(1);
        }
        gstream = &out_file;
    }
    for (int i = 0; i < ns->n; i++)
    {

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

            *gstream << std::setw(15) << std::right << ns->h_t[i];
            *gstream << std::setw(15) << std::right << ns->h_dt[i];
            *gstream << std::endl;
    }

    if(!print_screen)
    {
        out_file.close();
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

    if(print_screen)
    {
        gstream = &std::cout;
    }
    else
    {
        std::string ofname_log = ofname + ".log";
        out_file.open(ofname_log.c_str(), std::ios::out | std::ios::app );
        if(!out_file)
        {
          std::cerr << "gravidy: cannot open file "
                    << ofname_log.c_str()
                    << ": No such file or directory"
                    << std::endl;
          exit(1);
        }
        gstream = &out_file;
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

    if(!print_screen)
    {
        out_file.close();
    }
}
