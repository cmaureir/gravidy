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
#ifndef OPTIONSPARSER_HPP
#define OPTIONSPARSER_HPP
#include "../common.hpp"
#include <boost/program_options.hpp>
#include <sys/stat.h>

namespace po = boost::program_options;

/**
 * Class in charge to handle all the command-line parameters.
 * Besides handling the information, this class uses the boost/program_options
 * module to verify all the options.
 */

class OptionsParser {
    public:
        OptionsParser(int argc, char *argv[]);
        ~OptionsParser();

        std::string input_filename;
        std::string output_filename;
        std::string resume_filename;
        std::string snapshot_filename;
        float integration_time;
        float interval_time;
        int gpus;
        float snapshot_time;
        float softening;
        float eta;
        unsigned int snapshot_number;
        unsigned int length_output_number;
        bool resume;
        options ops;
        po::variables_map vm;
        po::options_description desc;

        bool file_exists(std::string filename);
        bool check_options();

    private:
        int argc;
        char *argv[];

};

#endif
