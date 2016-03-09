#include "Logger.hpp"

Logger::Logger(NbodySystem *ns)
{
    this->ns = ns;
    this->print_screen = this->ns->ops.print_screen;
    //this->ofname = this->get_timestamp() +"_"+ ns->output_filename;
    this->ofname = ns->output_filename;
    if (ns->resume)
    {
        this->ofname_info = ns->resume_filename;
    }
    else
    {
      this->ofname_info = ofname + ".info";
    }
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

void Logger::write_snapshot(int snapshot_number, double ITIME)
{
    std::ostringstream s;
    s << std::setw(4) << std::setfill('0') << snapshot_number;
    std::string ofname_all = ofname + ".snapshot_" + s.str();
    out_file.open(ofname_all.c_str(), std::ios::out);

    out_file << std::fixed;
    out_file.precision(4);
    out_file << "# Time:" <<std::setw(6) << ITIME;
    out_file << std::endl;

    for (int i = 0; i < ns->n; i++)
    {
        gstream->precision(2);
        out_file << std::setw(6)  << std::right << ns->h_id[i];

        // Scientific notation
        out_file << std::scientific;
        gstream->precision(6);
        out_file << std::setw(15) << std::right << ns->h_r[i].w;

        out_file << std::setw(15) << std::right << ns->h_r[i].x;
        out_file << std::setw(15) << std::right << ns->h_r[i].y;
        out_file << std::setw(15) << std::right << ns->h_r[i].z;

        out_file << std::setw(15) << std::right << ns->h_v[i].x;
        out_file << std::setw(15) << std::right << ns->h_v[i].y;
        out_file << std::setw(15) << std::right << ns->h_v[i].z;

        out_file << std::endl;
    }

    out_file.close();
}

void Logger::add_info(std::string key, std::string value)
{
    info_file.open(ofname_info.c_str(), std::ios_base::app);

    info_file << std::setw(2) << std::left  << "#";
    info_file << std::setw(20) << std::left  << key;
    info_file << std::setw(20) << std::right << value;
    info_file << std::endl;

    info_file.close();

}

bool file_exists(std::string filename)
{
    struct stat buffer;
    if (stat(filename.c_str(), &buffer) != -1)
    {
        return true;
    }
    return false;
}

void Logger::write_info()
{

    if(file_exists(ofname_info))
    {
        size_t fsize = ofname_info.size();
        std::string last = ofname_info.substr(fsize-1, fsize);
        bool is_number;
        int ext = 0;

        std::string newname = ofname_info + ".old";

        while(file_exists(newname))
        {
            // Check
            try
            {
                is_number = isdigit(last[0]);
            }
            catch(int e)
            {
                is_number = false;
            }

            if (is_number)
            {
                ext = (int)std::atoi(last.c_str()) + 1;
            }
            //newname =  ofname_info + ".old" + std::to_string(ext);// C++11
            newname =  ofname_info + ".old" + std::string(SSTR(ext));
            ext++;
        }

        int result = std::rename(ofname_info.c_str(), newname.c_str());
        if (result)
        {
            std::cerr << "gravidy: Failed to rename info file from "
                      << ofname_info
                      << " to "
                      << newname
                      << std::endl;
        }
    }
    info_file.open(ofname_info.c_str(), std::ios::out);


    info_file << std::setw(2) << std::left  << "#";
    info_file << std::setw(20) << std::left  << "NumberParticles:";
    info_file << std::setw(20) << std::right << ns->n;
    info_file << std::endl;

    info_file << std::setw(2) << std::left  << "#";
    info_file << std::setw(20) << std::left  << "Softening:";
    info_file << std::setw(20) << std::right << sqrt(ns->e2);
    info_file << std::endl;

    info_file << std::setw(2) << std::left  << "#";
    info_file << std::setw(20) << std::left  << "EtaTimestep:";
    info_file << std::setw(20) << std::right << ns->eta;
    info_file << std::endl;

    info_file << std::setw(2) << std::left  << "#";
    info_file << std::setw(20) << std::left  << "IntegrationTime:";
    info_file << std::setw(20) << std::right << ns->integration_time;
    info_file << std::endl;

    info_file << std::setw(2) << std::left  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "PrintScreen:";
    info_file << std::setw(20) << std::right << print_screen;
    info_file << std::endl;

    info_file << std::setw(2) << std::left  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "InputFilename:";
    info_file << std::setw(20) << std::right << ns->input_filename;
    info_file << std::endl;

    info_file << std::setw(2) << std::left  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "OutputFilename:";
    info_file << std::setw(20) << std::left << ns->output_filename;
    info_file << std::endl;

    info_file.close();
}

void Logger::print_info()
{

    if(print_screen)
    {
        gstream = &std::cout;
    }
    else
    {
        //std::string ofname_info = ofname + ".info";
        std::string ofname_info = ofname + ".log";
        out_file.open(ofname_info.c_str(), std::ios::out);
        gstream = &out_file;
    }

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
    *gstream << std::setw(8)  << std::right << ns->t_rlx;
    *gstream << std::endl;

    *gstream << std::setw(2)  << std::left  << "#";
    *gstream << std::setw(8)  << std::left  << "T_cr:";
    *gstream << std::setw(8)  << std::right << ns->t_cr;
    *gstream << std::endl;

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
        gstream = &out_file;
    }


    *gstream << "01 ";
    *gstream << std::fixed;
    *gstream << ITIME << " ";
    gstream->precision(4);
    for (int i = 0; i < 1/RADIUS_RATIO; i++) {
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
