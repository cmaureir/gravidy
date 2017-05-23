/*
 * Copyright (c) 2016
 *
 * Cristián Maureira-Fredes <cmaureirafredes@gmail.com>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "Logger.hpp"

Logger::Logger(NbodySystem *ns)
{
    this->ns = ns;
    this->print_screen = this->ns->ops.print_screen;
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

    s << "[";
    s << std::setw(4) << std::setfill('0') << 1900 + ltm->tm_year << "-";
    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_mon<< "-";
    s << std::setw(2) << std::setfill('0') << ltm->tm_mday << " ";

    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_hour << ":";
    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_min  << ":";
    s << std::setw(2) << std::setfill('0') << 1 + ltm->tm_sec  << "";
    s << "]";

    return s.str();
}

void Logger::write_snapshot(unsigned int snapshot_number, double ITIME)
{
    std::ostringstream s;
    s << std::setw(ns->length_snap) << std::setfill('0') << snapshot_number;
    std::string ofname_all = ofname + ".snapshot_" + s.str();
    out_file.open(ofname_all.c_str(), std::ios::out);

    out_file << std::fixed;
    out_file.precision(4);
    out_file << "# Time:" <<std::setw(6) << ITIME;
    out_file << std::endl;

    for (unsigned int i = 0; i < ns->n; i++)
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

    info_file << std::setw(2)  << std::left  << "#";
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


    info_file << std::setw(2)  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "NumberParticles:";
    info_file << std::setw(20) << std::right << ns->n;
    info_file << std::endl;

    info_file << std::setw(2)  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "Softening:";
    info_file << std::setw(20) << std::right << sqrt(ns->e2);
    info_file << std::endl;

    info_file << std::setw(2)  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "EtaTimestep:";
    info_file << std::setw(20) << std::right << ns->eta;
    info_file << std::endl;

    info_file << std::setw(2)  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "IntegrationTime:";
    info_file << std::setw(20) << std::right << ns->integration_time;
    info_file << std::endl;

    info_file << std::setw(2)  << std::left  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "PrintScreen:";
    info_file << std::setw(20) << std::right << print_screen;
    info_file << std::endl;

    info_file << std::setw(2)  << std::left  << std::left  << "#";
    info_file << std::setw(20) << std::left  << "InputFilename:";
    info_file << std::setw(20) << std::right << ns->input_filename;
    info_file << std::endl;

    info_file << std::setw(2)  << std::left  << std::left  << "#";
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
        std::string ofname_info = ofname + ".log";
        out_file.open(ofname_info.c_str(), std::ios::out);
        gstream = &out_file;
    }

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
    gstream->precision(6);

    size_t lr_size = sizeof(LAGRANGE_RADII)/4.0; // float in bytes is 4.0
    for (unsigned int i = 0; i < lr_size; i++)
    {
        *gstream << std::setw(6) << std::right << lagrange_radii[i] << " ";
    }
    *gstream << std::endl;

    if(!print_screen)
    {
        out_file.close();
    }
}

void Logger::print_all(double ITIME, unsigned int snapshot_number)
{
    if(print_screen)
    {
        gstream = &std::cout;
    }
    else
    {
        std::ostringstream s;
        s << std::setw(4) << std::setfill('0') << snapshot_number;
        std::string ofname_all = ofname + ".all.snap_" + s.str();
        out_file.open(ofname_all.c_str(), std::ios::out);
        gstream = &out_file;
    }
    for (unsigned int i = 0; i < ns->n; i++)
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

void Logger::print_energy_log(double ITIME, unsigned int iterations, long long interactions, unsigned int nsteps, double new_energy)
{
    ns->en.end = new_energy;
    double rel_error     = std::abs((ns->en.end - ns->en.tmp)/ns->en.tmp);
    double cum_error     = std::abs((ns->en.end - ns->en.ini)/ns->en.ini);
    //float grav_gflops = 0.0;
    ns->en.tmp = ns->en.end;

    float time = omp_get_wtime() - ns->gtime.integration_ini;
    if (ITIME != 0.0){
        ns->gtime.gflops =  KERNEL_GFLOP * (interactions / ns->gtime.update_end);
        //grav_gflops =  KERNEL_GFLOP * (interactions / ns->gtime.grav_end);
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
            //*gstream << std::setw(12)  << std::right << "GFLOPS(g)";
            *gstream << std::endl;
    }

    *gstream << std::fixed;
    *gstream << std::setw(2)  << std::right << "00";
    gstream->precision(3);
    *gstream << std::setw(10) << std::right << ITIME;
    *gstream << std::setw(15) << std::right << iterations;
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
    //*gstream << std::setw(12) << std::right << grav_gflops;
    *gstream << std::endl;

    if(!print_screen)
    {
        out_file.close();
    }
}

void Logger::log(int type,  std::string msg)
{
    // Error
    if (COLOR_OUTPUT)
    {
        if (type == 0)
        {
            std::cerr << color_blue << get_timestamp() << color_red << " [ERROR] " << color_disable << msg << std::endl;
        }
        // Warning
        else if (type == 1)
        {
            std::cerr << color_blue << get_timestamp() << color_yellow << " [WARNING] " << color_disable << msg << std::endl;
        }
        // Success
        else if (type == 2)
        {
            std::cerr << color_blue << get_timestamp() << color_green <<" [SUCCESS] " << color_disable << msg << std::endl;
        }
        // Information
        else if (type == 3)
        {
            std::cerr << color_blue << get_timestamp() << color_cyan << " [INFO] " << color_disable << msg << std::endl;
        }
    }
    else
    {
        if (type == 0)
        {
            std::cerr << get_timestamp() << " [ERROR] " << msg << std::endl;
        }
        // Warning
        else if (type == 1)
        {
            std::cerr << get_timestamp() << " [WARNING] " << msg << std::endl;
        }
        // Success
        else if (type == 2)
        {
            std::cerr << get_timestamp() << " [SUCCESS] " << msg << std::endl;
        }
        // Information
        else if (type == 3)
        {
            std::cerr << get_timestamp() << " [INFO] " << msg << std::endl;
        }
    }
}

void Logger::log_error(std::string msg)
{
    log(0, msg);
}

void Logger::log_warning(std::string msg)
{
    log(1, msg);
}

void Logger::log_success(std::string msg)
{
    log(2, msg);
}

void Logger::log_info(std::string msg)
{
    log(3, msg);
}
