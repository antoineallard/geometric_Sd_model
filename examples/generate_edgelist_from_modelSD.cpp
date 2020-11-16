/*
 *
 *  This command line code generates a random edgelist according to the S^D model.
 *
 *  Compilation requires at least the C++11 standard and must be done from the root repository of the project.
 *    Example: g++ -O3 -std=c++11 examples/generate_edgelist_from_modelSD.cpp -o generate_edgelist_from_modelSD
 *
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

#include "../src/modelSD_unix.hpp"

int main(int argc , char *argv[])
{
  // Initialize graph object.
  modelSD_t the_graph;

  // Parses the options and continues if everything is in order.
  if(parse_options(argc, argv, the_graph))
  {
    // Loads the hidden variables.
    the_graph.load_hidden_variables();

    // Generates an edgelist.
    the_graph.generate_graph();

    // Saves the graph and some metadada as a graphML file.
    the_graph.save_as_graphml();
  }

  // Returns successfully.
  return EXIT_SUCCESS;
}
