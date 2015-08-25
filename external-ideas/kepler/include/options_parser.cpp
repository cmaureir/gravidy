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
 * @fn check_options_noboost
 *
 * @param int argc
 *  Number of standard parameters
 * @param char* argv
 *  Parameters values
 *
 * @brief
 *   Checking the parameters of the program without using the
 *   boost library.
 *
 * @return state of the parameters
 */
bool check_options_noboost(int argc, char *argv[])
{
    input_file = std::string(argv[1]);
    itime = atoi(argv[2]);
    float e = atof(argv[3]);
    e2 = e*e;
    eta = atof(argv[4]);

    // Preparing output filename
    std::ostringstream ss;
    ss << "_t";
    ss << itime;
    ss << "_s";
    ss << e;
    ss << "_e";
    ss << eta;
    ss << ".out";
    std::string ext(ss.str());

    output_file = input_file+ext;

    return true;
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
    float e = 0.0;
    po::options_description desc("Options");
    desc.add_options()
        ("help,h",      "Display message")
        ("input,i",     po::value<std::string>(), "Input data filename")
        ("output,o",    po::value<std::string>(), "Output data filename")
        ("time,t",      po::value<float>(),       "Integration time")
        ("softening,s", po::value<float>(),       "Softening")
        ("eta,e",       po::value<float>(),       "ETA of time-step calculation")
        ("screen,p",    "Print summary in the screen instead of a file")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    /*
     * Help option
     */
    if (vm.count("help"))
    {
        std::cerr << desc << std::endl;
        return false;
    }

    /*
     * Input option
     */
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

    /*
     * Screen
     */
    if(vm.count("screen"))
    {
        print_log=0;
    }
    else
    {
        print_log = 1;
    }

    /*
     * Output option
     */
    std::ostringstream ss;
    ss << "_t";
    ss << itime;
    ss << "_s";
    ss << e;
    ss << "_e";
    ss << eta;
    ss << ".out";
    std::string ext(ss.str());
    if (vm.count("output"))
    {
        output_file = vm["output"].as<std::string>();
        output_file = output_file+ext;
    }
    else
    {
        output_file = input_file+ext;
    }

    /*
     * Time option
     */
    if (vm.count("time"))
    {
        itime = vm["time"].as<float>();
    }
    else
    {
        std::cerr << "gravidy: option requires an argument -- 'time'"
                  << std::endl;
        std::cerr << desc
                  << std::endl;
        return false;
    }

    /*
     * Softening
     */
    e2 = E*E;
    if (vm.count("softening"))
    {
        e = vm["softening"].as<float>();
        e2 = e*e;
    }

    /*
     * Eta
     */
    eta  = ETA_N;
    if (vm.count("eta"))
    {
        eta = vm["eta"].as<float>();
    }


    return true;
}
