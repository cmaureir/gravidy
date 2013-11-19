#ifndef OPTIONSPARSER_HPP
#define OPTIONSPARSER_HPP
#include "../common.hpp"
#include <boost/program_options.hpp>
#include <sys/stat.h>
#include <string>

namespace po = boost::program_options;

class OptionsParser {
    public:
        OptionsParser(int argc, char *argv[]);
        ~OptionsParser();

        std::string input_filename;
        std::string output_filename;
        float integration_time;
        float softening;
        float eta;
        int print_log;
        po::variables_map vm;
        po::options_description desc;

        bool file_exists(std::string filename);
        bool check_options();

    private:
        int argc;
        char *argv[];

};

#endif
