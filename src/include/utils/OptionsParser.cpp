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
#include "OptionsParser.hpp"

/** Constructor which includes the options and requirements for the
 * command line interface
 */
OptionsParser::OptionsParser(int argc, char *argv[])
{

    po::options_description help("GraviDy");
    help.add_options()("help,h", "Display this message");

    // Required
    po::options_description main("Required options");
    main.add_options()
        ("input,i", po::value<std::string>()->value_name("<filename>"),
           "Input data filename")
        ("time,t", po::value<float>()->value_name("<value>")->default_value(1),
            "Integration time (In N-body units)")
        ("resume,r", po::value<std::string>()->value_name("<filename>.info"),
           "Resume a simulation with an .info file")
    ;

    // Elective
    po::options_description elective("Elective options");
    elective.add_options()
        ("output,o", po::value<std::string>()->value_name("<filename>"),
            "Output data filename")
        ("softening,s", po::value<float>()->value_name("<value>"),
            "Softening parameter (default 1e-4)")
        ("eta,e", po::value<float>()->value_name("<value>"),
            "ETA of time-step calculation (default 0.01)")
        ("screen,p", "Print summary in the screen instead of a file")
    ;

    // Extra
    po::options_description extra("Extra options");
    extra.add_options()
        ("lagrange,l",
            "Print information of the Lagrange Radii in every integration time")
        ("all,a",
            "Print all the information of N-particles in every integration time")
        ("gpu,g", po::value<int>()->value_name("<value>")->default_value(0),
            "GPUs to use, by default is the maximum available devices (use even numbers)")
    ;

    desc.add(help);
    desc.add(main);
    desc.add(elective);
    desc.add(extra);
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
}

/** Destructor */
OptionsParser::~OptionsParser()
{
    // Empty destructor
}

/** Method that check if a file exist on the system */
bool OptionsParser::file_exists(std::string filename)
{
    struct stat buffer;
    if (stat(filename.c_str(), &buffer) != -1)
    {
        return true;
    }
    return false;
}

/** Method in charge to verify all the command-line arguments that where
 * used for the current run
 */
bool OptionsParser::check_options()
{

    if (vm.count("help"))
    {
        std::cerr << desc << std::endl;
        return false;
    }

    double config_time = 0.0;
    snapshot_time = 0.0;
    if (vm.count("resume"))
    {

        resume_filename = vm["resume"].as<std::string>();
        resume = true;


        if(!file_exists(resume_filename))
        {
            std::cerr << "gravidy: cannot access "
                      << resume_filename
                      << ": No such file or directory"
                      << std::endl;
            return false;
        }
        else
        {
            std::cerr << "Reading resume_filename: " << resume_filename << std::endl;

            // Open resume_filename
            std::ifstream rfile(resume_filename.c_str());
            std::string key;
            std::string value;
            std::map<std::string,std::string> config_param;
            std::set<std::string> allowed_param;

            allowed_param.insert("NumberParticles");
            allowed_param.insert("Softening");
            allowed_param.insert("EtaTimestep");
            allowed_param.insert("IntegrationTime");
            allowed_param.insert("PrintScreen");
            allowed_param.insert("InputFilename");
            allowed_param.insert("OutputFilename");
            allowed_param.insert("SnapshotNumber");

            if (rfile.is_open())
            {
                while (rfile.good())
                {

                    // Check EOF
                    if(rfile.eof())
                        break;

                    // Get key
                    if (getline(rfile, key, ':'))
                    {
                        // Clean key
                        char kchars[] = "# ";
                        for (unsigned int i = 0; i < strlen(kchars); i++)
                        {
                            key.erase(remove(key.begin(), key.end(), kchars[i]),
                                key.end());
                        }

                        if (getline(rfile, value))
                        {
                            // Clean value
                            char vchars[] = " ";
                            for (unsigned int i = 0; i < strlen(vchars); i++)
                            {
                                value.erase(remove(value.begin(), value.end(),
                                    vchars[i]), value.end());
                            }
                            config_param[key] = value;
                        }
                    }
                }
            }
            // Closing resume file
            rfile.close();

            // Check is all the `key` are valid
            std::map<std::string, std::string>::iterator ii;
            for(ii=config_param.begin(); ii != config_param.end(); ii++)
            {
                if (allowed_param.find((*ii).first) != allowed_param.end())
                {
                    std::cerr << (*ii).first << " -> " << (*ii).second << std::endl;
                }
                else
                {
                    std::cerr << "Invalid param" << std::endl;
                    std::cerr << "gravidy: Invalid option -- "
                              << (*ii).first
                              << " in the info file"
                              << std::endl;
                    std::cerr << desc
                              << std::endl;
                    return false;
                }
            }

            // Check if all the `value` of the `key` are valid

            // Checking time
            if (vm.count("time"))
            {
                integration_time = vm["time"].as<float>();
                config_time = strtod(config_param["IntegrationTime"].c_str(), NULL);
                if (config_time >= integration_time)
                {
                    std::cerr << "gravidy: option requires an argument -- 'time'"
                              << " greater than the one specifid in the info file"
                              << std::endl;
                    std::cerr << desc
                              << std::endl;
                    return false;
                }
                else
                {
                    snapshot_time = config_time;
                }
            }
            else
            {
                std::cerr << "gravidy: option requires an argument -- 'time'"
                          << std::endl;
                std::cerr << desc
                          << std::endl;
                return false;
            }

            // Checking softening
            softening = strtod(config_param["Softening"].c_str(), NULL);
            if (vm.count("softening"))
            {
                softening = vm["softening"].as<float>();
            }

            // Checking eta
            eta = strtod(config_param["EtaTimestep"].c_str(), NULL);
            if (vm.count("eta"))
            {
                eta = vm["eta"].as<float>();
            }

            // Checking PrintScreen
            //ops.print_screen = std::stoi(config_param["PrintScreen"].c_str()); //C++11
            ops.print_screen = std::atoi(config_param["PrintScreen"].c_str());
            if (vm.count("PrintScreen"))
            {
                ops.print_screen=1;
            }

            // check NumberParticles ?
            //  ...Check the number of lines
            //  ...TODO what if some particles are removed in the future?

            // check InputFilename
            // keep the same structure of the files.
            std::string InputFilename = config_param["InputFilename"];
            input_filename = InputFilename;
            std::ostringstream ss;
            ss << ".out";
            std::string ext(ss.str());
            output_filename = InputFilename+ext;

            // check OutputFilename
            std::string OutputFilename = config_param["OutputFilename"];

            // Checking snapshot number
            //snapshot_number = std::stoi(config_param["SnapshotNumber"].c_str()); // C++11
            snapshot_number = std::atoi(config_param["SnapshotNumber"].c_str());

            std::ostringstream s;
            s << std::setw(4) << std::setfill('0') << snapshot_number;
            std::string snap_name = OutputFilename + ".snapshot_" + s.str();

            if (!file_exists(snap_name))
            {
                std::cerr << "gravidy: cannot access "
                          << snap_name
                          << ": No such file or directory\n"
                          << "(Check the OutputFilename and SnapshotNumber)"
                          << std::endl;
                return false;
            }
            else
            {
                snapshot_filename = snap_name;
                //input_filename = snap_name;
            }
        }
    }
    else // if the simulation is not being resume
    {
        resume = false;
        snapshot_number = 0;
        if (vm.count("input"))
        {
            input_filename = vm["input"].as<std::string>();

            if(!file_exists(input_filename))
            {
                std::cerr << "gravidy: cannot access "
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

        if(vm.count("gpu"))
            gpus = vm["gpu"].as<int>();

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


        if (vm.count("output"))
        {
            std::ostringstream ss;
            ss << ".out";
            std::string ext(ss.str());
            output_filename = vm["output"].as<std::string>();
            output_filename = output_filename+ext;
        }
        else
        {

            std::ostringstream ss;
            ss << ".out";
            std::string ext(ss.str());
            output_filename = input_filename+ext;
        }
    }


    return true;
}
