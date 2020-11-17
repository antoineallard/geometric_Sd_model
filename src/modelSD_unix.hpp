/*
 *
 *  Provides the functions related to UNIX operating system.
 *
 *  Author:  Antoine Allard
 *  WWW:     antoineallard.info
 *  Date:    November 2020
 *
 *  Copyright (C) 2020  Antoine Allard
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef MODELSD_UNIX_HPP_INCLUDED
#define MODELSD_UNIX_HPP_INCLUDED

// Standard Template Library
#include <cstdlib>
#include <string>
// Operating System
#include <unistd.h>
// embeddingS1
#include "modelSD.hpp"





// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Prints the information on the way the command line UNIX program should be used.
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void print_usage()
{
  std::cout << std::endl;
  std::cout << "NAME"                                                                                 << std::endl;
  std::cout << "\tgeneratingSD -- a program to generate complex networks in the SD metric space"      << std::endl;
  std::cout << "SYNOPSIS"                                                                             << std::endl;
  std::cout << "\tmodelSD [options] <hidden_variables_filename>"                                 << std::endl;
  std::cout << "DESCRIPTION"                                                                          << std::endl;
  std::cout << "\tSee README.md for further details."                                                 << std::endl;
  std::cout << std::endl;
}





// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Prints the information about the options of the command line UNIX program should be used.
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void print_help()
{
  std::cout << std::endl;
  std::cout << "The following options are available:"                                                                                                                     << std::endl;
  std::cout << "\t-b [VALUE]     Specifies the value for parameter BETA."                                                                                                 << std::endl;
  std::cout << "\t-d [VALUE]     Specifies the DIMENSION of the hypersphere. Default: DIMENSION = 1."                                                                     << std::endl;
  std::cout << "\t-h             Print this message on screen and exit."                                                                                                  << std::endl;
  std::cout << "\t-m [VALUE]     Specifies the value for parameter MU."                                                                                                   << std::endl;
  std::cout << "\t               Default: MU = BETA * sin(DIMENSION*PI/BETA) * gamma(DIMENSION/2) / ( 2 * PI * average_kappa * pow(PI, DIMENSION/2) )."                   << std::endl;
  std::cout << "\t-n             Indicates that the first column of the hidden variables file provides the name of the vertices. Otherwise, vertices are"                 << std::endl;
  std::cout << "\t               assigned a name following the convention: v####."                                                                                        << std::endl;
  std::cout << "\t-o [ROOTNAME]  Specifies the rootname used for all output files. Uses the filename of the hidden variables file as rootname if not specified."          << std::endl;
  std::cout << "\t-s [SEED]      Program uses a custom seed for the random number generator. Default: EPOCH."                                                             << std::endl;
  std::cout << "\t-t             Indicates that the last columns of the hidden variables file provides the angular positions of the vertices."                            << std::endl;
  std::cout << "\t-v             Outputs the hidden variables (kappa and angular position) used to the generate the network into a file (uses the edgelist's rootname)."  << std::endl;
  std::cout << std::endl;
}





// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Parses the options (for UNIX-like command line use) and returns the filename of the edgelist or
//   a flag indicating that help dialogues have been shown on screen and further computation is not
//   required.
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
bool parse_options(int argc , char *argv[], pgl::modelSD_t &the_graph)
{
  // Shows the options if no argument is given.
  if(argc == 1)
  {
    print_usage();
    print_help();
    return false;
  }

  // <hidden_variables_filename>
  the_graph.HIDDEN_VARIABLES_FILENAME = argv[argc - 1];

  // Parsing options.
  int opt;
  while ((opt = getopt(argc,argv,"b:d:hm:no:s:tv")) != -1)
  {
    switch(opt)
    {
      case 'b':
        the_graph.BETA = std::stod(optarg);
        break;

      case 'd':
        the_graph.DIMENSION = std::stod(optarg);
        break;

      case 'h':
        print_usage();
        print_help();
        return false;

      case 'm':
        the_graph.MU = std::stod(optarg);
        break;

      case 'n':
        the_graph.NAME_PROVIDED = true;
        break;

      case 'o':
        the_graph.OUTPUT_ROOTNAME = optarg;
        the_graph.CUSTOM_OUTPUT_ROOTNAME_MODE = true;
        break;

      case 's':
        the_graph.SEED = std::stoi(optarg);
        the_graph.initialize_random_number_generator();
        break;

      case 't':
        the_graph.THETA_PROVIDED = true;
        break;

      case 'v':
        the_graph.OUTPUT_VERTICES_PROPERTIES = true;
        break;

      default:
        print_usage();
        print_help();
        return false;
    }
  }

  // Uses the default rootname for output files.
  if(the_graph.CUSTOM_OUTPUT_ROOTNAME_MODE == false)
  {
    size_t lastdot = the_graph.HIDDEN_VARIABLES_FILENAME.find_last_of(".");
    if(lastdot == std::string::npos)
    {
      the_graph.OUTPUT_ROOTNAME = the_graph.HIDDEN_VARIABLES_FILENAME;
    }
    the_graph.OUTPUT_ROOTNAME = the_graph.HIDDEN_VARIABLES_FILENAME.substr(0, lastdot);
  }

  // Indicates that everything is in order.
  return true;
}


#endif // MODELSD_UNIX_HPP_INCLUDED
