#include "OptionsParser.hpp"

OptionsParser::OptionsParser(int argc, char *argv[])
{
    desc.add_options()
        ("help,h",      "Display message")
        ("input,i",     po::value<std::string>(), "Input data filename")
        ("output,o",    po::value<std::string>(), "Output data filename")
        ("time,t",      po::value<float>(),       "Integration time")
        ("softening,s", po::value<float>(),       "Softening")
        ("eta,e",       po::value<float>(),       "ETA of time-step calculation")
        ("screen,p",    "Print summary in the screen instead of a file")
    ;

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
}

OptionsParser::~OptionsParser()
{
}

bool OptionsParser::file_exists(std::string filename)
{
    struct stat buffer;
    if (stat(filename.c_str(), &buffer) != -1)
    {
        return true;
    }
    return false;
}


bool OptionsParser::check_options()
{

    if (vm.count("help"))
    {
        std::cerr << desc << std::endl;
        return false;
    }

    if (vm.count("input"))
    {
        input_filename = vm["input"].as<std::string>();

        if(!file_exists(input_filename))
        {
            std::cout << "gravidy: cannot access "
                      << input_filename
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

    if(vm.count("screen"))
        print_log=0;
    else
        print_log = 1;


    if (vm.count("time"))
        integration_time = vm["time"].as<float>();
    else
    {
        std::cerr << "gravidy: option requires an argument -- 'time'"
                  << std::endl;
        std::cerr << desc
                  << std::endl;
        return false;
    }

    softening = E;
    if (vm.count("softening"))
    {
        softening = vm["softening"].as<float>();
    }

    eta  = ETA_N;
    if (vm.count("eta"))
        eta = vm["eta"].as<float>();

    std::ostringstream ss;
    ss << "_t";
    ss << integration_time;
    ss << "_s";
    ss << softening;
    ss << "_e";
    ss << eta;
    ss << ".out";
    std::string ext(ss.str());

    if (vm.count("output"))
    {
        output_filename = vm["output"].as<std::string>();
        output_filename = output_filename+ext;
    }
    else
    {
        output_filename = input_filename+ext;
    }


    return true;
}
