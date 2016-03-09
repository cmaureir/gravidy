#ifndef LOGGER_HPP
#define LOGGER_HPP
#include "../common.hpp"
#include "../NbodySystem.hpp"
#include <cctype>
#include <sstream>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

class Logger {
    public:
        Logger(NbodySystem *ns);
        ~Logger();

        NbodySystem *ns;

        std::ofstream out_file;
        std::ofstream info_file;
        std::ostream *gstream;
        bool print_screen;
        std::string ofname;
        std::string ofname_info;

        std::string get_timestamp();
        void print_info();
        void add_info(std::string key, std::string value);
        void write_info();
        void write_snapshot(int snapshot_number, double ITIME);
        void print_all(double ITIME);
        void print_energy_log(double ITIME, int iterations, long long interactions,
                              int nsteps, double new_energy);
        void print_lagrange_radii(double ITIME, std::vector<double> lagrange_radii);
};

#endif // LOGGER_HPP
