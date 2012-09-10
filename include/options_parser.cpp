#include "options_parser.hpp"

/*
 * @fn file_exist(const std::string)
 *
 * @brief
 *  Obtained from
 *  http://stackoverflow.com/questions/4316442/stdofstream-check-if-file-exists-before-writing
 * @param filename
 *  File name to check
 *
 * @return True if exist, else False
 */
bool file_exists(const std::string& filename)
{
    struct stat buffer;
    if (stat(filename.c_str(), &buffer) != -1)
    {
        return true;
    }
    return false;
}

/*
 * @fn check_options
 *
 * @param int argc
 *  Number of standard parameters
 * @param char* argv
 *  Parameters values
 *
 * @brief
 *   Checking the parameters of the program using the
 *   boost library.
 *
 * @return state of the parameters
 */
bool check_options(int argc, char *argv[])
{
    po::options_description desc("Options");
    desc.add_options()
        ("help,h"   , "Display message")
        ("input,i"  , po::value<std::string>() , "Input data filename")
        ("output,o" , po::value<std::string>() , "Output data filename")
        ("time,t"   , po::value<double>()      , "Integration time")
        ("run,r"    , po::value<std::string>()  , "Running type, CPU (default) or GPU")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Help option
    if (vm.count("help"))
    {
        std::cerr << desc
                  << std::endl;
        return false;
    }

    // Input option
    if (vm.count("input"))
    {
        input_file = vm["input"].as<std::string>();
        if(!file_exists(input_file))
        {
            std::cout << "gravidy: cannot access "
                      << input_file
                      << ": No such file or directory"
                      << std::endl;
            return false;
        }
    }
    else
    {
        std::cerr << "gravidy: option requires an argument -- 'input'"
                  << std::endl;
        std::cerr << desc
                  << std::endl;
        return false;
    }

    // Output option
    if (vm.count("output"))
    {
        output_file = vm["output"].as<std::string>();
    }
    //else
    //{
    //    std::cerr << "gravidy: option requires an argument -- 'output'" << std::endl;
    //    std::cerr << desc << std::endl;
    //    return false;
    //}

    // Time option
    if (vm.count("time"))
    {
        int_time = vm["time"].as<double>();
    }
    else
    {
        std::cerr << "gravidy: option requires an argument -- 'time'"
                  << std::endl;
        std::cerr << desc
                  << std::endl;
        return false;
    }

    // Run option
    run = std::string("cpu"); // Default
    if (vm.count("run"))
    {
        run = vm["run"].as<std::string>();
        if(!(run == "gpu" || run == "cpu"))
        {
            std::cerr << "gravidy: Invalid option -- '"
                      << run
                      << "'"
                      << std::endl;
            return false;
        }
    }

    return true;
}
