#include <iostream>
#include <boost/program_options.hpp>
#include <sys/stat.h>
#ifndef COMMON_HPP
#define COMMON_HPP
#include "common.hpp"
#endif

namespace po = boost::program_options;
bool check_options_noboost(int argc, char *argv[]);
bool check_options(int argc, char *argv[]);
