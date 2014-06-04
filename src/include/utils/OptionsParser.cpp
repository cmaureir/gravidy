#include "OptionsParser.hpp"

OptionsParser::OptionsParser(int argc, char *argv[])
{


    po::options_description help("GraviDy");
    help.add_options()("help,h",      "Display this message");

    po::options_description main("Required options");
    main.add_options()
        //("input,i",     po::value<std::string>()->value_name("<filename>"), "Input data filename")
        //("time,t",      po::value<float>()->value_name("<value>")->default_value(1),       "Integration time (In N-body units)")
        ("input,i",     po::value<std::string>(), "Input data filename")
        ("time,t",      po::value<float>(),       "Integration time (In N-body units)")
    ;

    po::options_description optional("Optional options");
    optional.add_options()
        //("output,o",    po::value<std::string>()->value_name("<filename>"), "Output data filename")
        //("softening,s", po::value<float>()->value_name("<value>"),       "Softening parameter (default 1e-4)")
        //("eta,e",       po::value<float>()->value_name("<value>"),       "ETA of time-step calculation (default 0.01)")
        ("output,o",    po::value<std::string>(), "Output data filename")
        ("softening,s", po::value<float>(),       "Softening parameter (default 1e-4)")
        ("eta,e",       po::value<float>(),       "ETA of time-step calculation (default 0.01)")
        ("screen,p",    "Print summary in the screen instead of a file")
    ;

    po::options_description extra("Extra options");
    extra.add_options()
        ("lagrange,l",    "Print information of the Lagrange Radii in every integration time")
        ("all,a",    "Print all the information of all the particles in every integration time")
    ;

    desc.add(help);
    desc.add(main);
    desc.add(optional);
    desc.add(extra);
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

    // Options structure
    ops.print_screen = 0;
    ops.print_all = 0;
    ops.print_lagrange = 0;

    if(vm.count("screen"))
        ops.print_screen=1;
    if(vm.count("all"))
        ops.print_all = 1;
    if(vm.count("lagrange"))
        ops.print_lagrange = 1;

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
