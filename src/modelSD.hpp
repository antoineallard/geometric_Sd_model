/*
 *
 *  This class provides the functions to generate a graph in the S^D space.
 *
 *  Compilation requires at least the C++11 standard.
 *    Example: g++ -O3 -std=c++11 my_code.cpp -o my_program
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

#ifndef MODELSD_HPP_INCLUDED
#define MODELSD_HPP_INCLUDED

// Standard Template Library
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>



class modelSD_t
{
  // Flags controlling options.
  public:
    bool CUSTOM_OUTPUT_ROOTNAME_MODE = false;
    bool NAME_PROVIDED = false;
    // bool NATIVE_INPUT_FILE = false;
    bool THETA_PROVIDED = false;
    bool OUTPUT_VERTICES_PROPERTIES = false;
  // Global parameters.
  public:
    // Random number generator seed.
    int SEED = std::time(NULL);
    // Parameter beta (clustering).
    double BETA = -1;
    // Number of dimensions.
    int DIMENSION = 1;
    // Parameter mu (average degree).
    double MU = -1;
    // Radius of DIMENSION-sphere.
    double RADIUS = 0;
    // Rootname for the output files;
    std::string OUTPUT_ROOTNAME = "default_output_rootname";
    // Input hidden variables filename.
    std::string HIDDEN_VARIABLES_FILENAME;
  // General internal objects.
  private:
    // pi
    const double PI = 3.141592653589793238462643383279502884197;
    // Random number generator
    std::mt19937 engine;
    std::uniform_real_distribution<double> uniform_01;
    // Mapping the numerical ID of vertices to their name.
    std::vector<std::string> Num2Name;
  // Objects related to the graph ensemble.
  private:
    // Number of vertices.
    int nb_vertices = 0;
    // Hidden variables of the vertices.
    std::vector<double> kappa;
    // Positions of the vertices.
    std::vector< std::vector<double> > theta;
    // Edgelist
    std::set< std::pair<int, int> > edgelist;
    // Vectors containing the expected and real degrees.
    std::vector<double> edegree;
    std::vector<int> rdegree;
  // Public functions to generate the graphs.
  public:
    // Constructor.
    modelSD_t() { initialize_random_number_generator(); };
    // Loads the values of the hidden variables (i.e., kappa and angular positions).
    void load_hidden_variables();
    // Generates an edgelist and writes it into a file.
    void generate_graph();
    // Random angular positions (if not provided).
    void generate_random_angular_positions();
    // Initializes the random number generator to the seed using SEED.
    void initialize_random_number_generator();
    // Saves the graph and (some metadada in another file) in an edgelist file.
    void save_as_edgelist(int width = 15);
    // Saves the graph and some metadada in a graphML file.
    void save_as_graphml();
    // Sets MU to its default value (obtained in the limit N -> infinity).
    void set_mu_to_default_value();
  // Private functions linked to the generation of a random edgelist.
  private:
    // Computes the radius of the hypersphere.
    void compute_radius();
    // Connection probability.
    double compute_connection_probability(int v1, int v2);
    // Saves the values of the hidden variables (i.e., kappa and theta).
    void save_vertices_properties(std::vector<int>& rdegree, std::vector<double>& edegree, int width);
    // Gets and format current date/time.
    std::string get_time();
};





// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
double modelSD_t::compute_connection_probability(int v1, int v2)
{
  double angular_distance = 0;
  if(DIMENSION == 1)
  {
    angular_distance = PI - std::fabs(PI - std::fabs(theta[v1][0] - theta[v2][0]));
    return 1 / (1 + std::pow((RADIUS * angular_distance) / (MU * kappa[v1] * kappa[v2]), BETA));
  }
  else if(DIMENSION == 2)
  {
    angular_distance = std::cos(std::fabs(theta[v1][0] - theta[v2][0]));
    angular_distance *= std::sin(theta[v1][1]) * std::sin(theta[v2][1]);
    angular_distance += std::cos(theta[v1][1]) * std::cos(theta[v2][1]);
    angular_distance = std::acos(angular_distance);
    return 1. / (1. + std::pow((RADIUS * angular_distance) / std::sqrt(MU * kappa[v1] * kappa[v2]), BETA));
  }
  else
  {
    // general formula
    // return 1 / (1 + std::pow((RADIUS * angular_distance) / std::pow(MU * kappa[v1] * kappa[v2], 1. / DIMENSION), BETA));
    std::cerr << "ERROR: No angular distance has been implemented for DIMENSION > 2." << std::endl;
    std::terminate();
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::compute_radius()
{
  double D = DIMENSION;
  RADIUS = nb_vertices * std::tgamma((D + 1) / 2);
  RADIUS /= 2 * std::pow(PI, (D + 1) / 2);
  RADIUS = std::pow(RADIUS, 1. / D);
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::generate_graph()
{
  // Initializes the containers.
  edgelist.clear();
  edegree.clear();
  rdegree.clear();
  // Initializes the containers for the expected and real degrees.
  if(OUTPUT_VERTICES_PROPERTIES)
  {
    edegree.resize(nb_vertices, 0);
    rdegree.resize(nb_vertices, 0);
  }
  // Makes sure the value of beta has been provided.
  if(BETA < 0)
  {
    std::cerr << "ERROR: The value of parameter beta must be provided." << std::endl;
    std::terminate();
  }
  // Sets the value of mu, if not provided.
  if(MU < 0)
  {
    set_mu_to_default_value();
  }
  // Generates the edgelist.
  double prob;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    for(int v2(v1 + 1); v2<nb_vertices; ++v2)
    {
      prob = compute_connection_probability(v1, v2);
      if(uniform_01(engine) < prob)
      {
        edgelist.insert(std::make_pair(v1, v2));
        if(OUTPUT_VERTICES_PROPERTIES)
        {
          rdegree[v1] += 1;
          rdegree[v2] += 1;
        }
      }
      if(OUTPUT_VERTICES_PROPERTIES)
      {
        edegree[v1] += prob;
        edegree[v2] += prob;
      }
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::generate_random_angular_positions()
{
  theta.clear();
  theta.resize(nb_vertices, std::vector<double>(DIMENSION));
  for(int v(0); v<nb_vertices; ++v)
  {
    if(DIMENSION == 1)
    {
      theta[v][0] = 2 * PI * uniform_01(engine);
    }
    else if(DIMENSION == 2)
    {
      theta[v][0] = 2 * PI * uniform_01(engine);
      theta[v][1] = std::acos(2 * uniform_01(engine) - 1);
    }
    else
    {
      std::cerr << "ERROR: Could not generate random angular positions for DIMENSION > 2." << std::endl;
      std::terminate();
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::string modelSD_t::get_time()
{
  // Gets the current date/time.
  time_t theTime = time(NULL);
  struct tm *aTime = gmtime(&theTime);
  int year    = aTime->tm_year + 1900;
  int month   = aTime->tm_mon + 1;
  int day     = aTime->tm_mday;
  int hours   = aTime->tm_hour;
  int minutes = aTime->tm_min;
  // Format the string.
  std::string the_time = std::to_string(year) + "/";
  if(month < 10)
    the_time += "0";
  the_time += std::to_string(month) + "/";
  if(day < 10)
    the_time += "0";
  the_time += std::to_string(day) + " " + std::to_string(hours) + ":";
  if(minutes < 10)
    the_time += "0";
  the_time += std::to_string(minutes) + " UTC";
  // Returns the date/time.
  return the_time;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::initialize_random_number_generator()
{
  // Initializes the random number generator.
  engine.seed(SEED);
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::load_hidden_variables()
{
  // Resets the number of vertices.
  nb_vertices = 0;
  // Resets the containers.
  kappa.clear();
  theta.clear();
  Num2Name.clear();
  // Stream object.
  std::stringstream one_line;
  // String objects.
  std::string full_line, name1_str, name2_str;
  // Opens the stream and terminates if the operation did not succeed.
  std::fstream hidden_variables_file(HIDDEN_VARIABLES_FILENAME.c_str(), std::fstream::in);
  if( !hidden_variables_file.is_open() )
  {
    std::cerr << "Could not open file: " << HIDDEN_VARIABLES_FILENAME << "." << std::endl;
    std::terminate();
  }
  // Reads the hidden variables file line by line.
  while( !hidden_variables_file.eof() )
  {
    // Reads the first entry of a line of the file.
    std::getline(hidden_variables_file, full_line);
    hidden_variables_file >> std::ws;
    one_line.str(full_line);
    one_line >> std::ws;
    one_line >> name1_str >> std::ws;
    // Skips lines of comment.
    if(name1_str == "#")
    {
      one_line.clear();
      continue;
    }
    // Adds the new vertex and its hidden variable(s).
    if(NAME_PROVIDED)
    {
      one_line >> name2_str >> std::ws;
      Num2Name.push_back(name1_str);
      kappa.push_back(std::stod(name2_str));
    }
    else
    {
      Num2Name.push_back("v" + std::to_string(nb_vertices));
      kappa.push_back(std::stod(name1_str));
    }
    if(THETA_PROVIDED)
    {
      theta.push_back(std::vector<double>(DIMENSION));
      for(int d(0); d<DIMENSION; ++d)
      {
        one_line >> name2_str >> std::ws;
        theta.back()[d] = std::stod(name2_str);
      }
    }
    ++nb_vertices;
    one_line.clear();
  }
  // Closes the stream.
  hidden_variables_file.close();
  // Generates the angular positions, if not provided.
  if(!THETA_PROVIDED)
  {
    generate_random_angular_positions();
  }
  // Computes the radius of the D-sphere.
  compute_radius();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::save_as_edgelist(int width)
{
  // Sets the name of the file to write the edgelist into.
  std::string edgelist_filename = OUTPUT_ROOTNAME + "_edgelist.dat";
  // Opens the stream and terminates if the operation did not succeed.
  std::fstream edgelist_file(edgelist_filename.c_str(), std::fstream::out);
  if( !edgelist_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << edgelist_filename << "." << std::endl;
    std::terminate();
  }
  // Writes the header.
  edgelist_file << "# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=" << std::endl;
  edgelist_file << "# Generated on:           " << get_time()                << std::endl;
  edgelist_file << "# Hidden variables file:  " << HIDDEN_VARIABLES_FILENAME << std::endl;
  edgelist_file << "# Seed:                   " << SEED                      << std::endl;
  edgelist_file << "#"                                                       << std::endl;
  edgelist_file << "# Parameters"                                            << std::endl;
  edgelist_file << "#   - nb. vertices:       " << nb_vertices               << std::endl;
  edgelist_file << "#   - dimension:          " << DIMENSION                 << std::endl;
  edgelist_file << "#   - beta:               " << BETA                      << std::endl;
  edgelist_file << "#   - mu:                 " << MU                        << std::endl;
  edgelist_file << "#   - radius:             " << RADIUS                    << std::endl;
  edgelist_file << "# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=" << std::endl;
  edgelist_file << "#";
  edgelist_file << std::setw(width - 1) << "Vertex1" << " ";
  edgelist_file << std::setw(width)     << "Vertex2" << " ";
  edgelist_file << std::endl;

  for(auto edge : edgelist)
  {
    edgelist_file << std::setw(width) << Num2Name[edge.first] << " ";
    edgelist_file << std::setw(width) << Num2Name[edge.second] << " ";
    edgelist_file << std::endl;
  }
  // Closes the stream.
  edgelist_file.close();
  // Outputs the hidden variables, if required.
  if(OUTPUT_VERTICES_PROPERTIES)
  {
    save_vertices_properties(rdegree, edegree, width);
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::save_as_graphml()
{
  // Sets the name of the file to write the edgelist into.
  std::string graphml_filename = OUTPUT_ROOTNAME + ".xml";
  // Opens the stream and terminates if the operation did not succeed.
  std::fstream graphml_file(graphml_filename.c_str(), std::fstream::out);
  if( !graphml_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << graphml_filename << "." << std::endl;
    std::terminate();
  }

  // Header
  graphml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"                                                                                   << std::endl;
  graphml_file << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\""                                                                     << std::endl;
  graphml_file << "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""                                                             << std::endl;
  graphml_file << "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">" << std::endl;

  // Declares the properties
  graphml_file << std::endl;
  graphml_file << "  <!-- property keys -->"                                                                      << std::endl;
  graphml_file << "  <key id=\"key0\" for=\"graph\" attr.name=\"generated_on\" attr.type=\"string\" />"           << std::endl;
  graphml_file << "  <key id=\"key1\" for=\"graph\" attr.name=\"hidden_variables_file\" attr.type=\"string\" />"  << std::endl;
  graphml_file << "  <key id=\"key2\" for=\"graph\" attr.name=\"seed\" attr.type=\"int\" />"                      << std::endl;
  graphml_file << "  <key id=\"key3\" for=\"graph\" attr.name=\"dimension\" attr.type=\"int\" />"                 << std::endl;
  graphml_file << "  <key id=\"key4\" for=\"graph\" attr.name=\"beta\" attr.type=\"double\" />"                   << std::endl;
  graphml_file << "  <key id=\"key5\" for=\"graph\" attr.name=\"mu\" attr.type=\"double\" />"                     << std::endl;
  graphml_file << "  <key id=\"key6\" for=\"graph\" attr.name=\"radius\" attr.type=\"double\" />"                 << std::endl;
  graphml_file << "  <key id=\"key7\" for=\"node\" attr.name=\"name\" attr.type=\"string\" />"                    << std::endl;
  graphml_file << "  <key id=\"key8\" for=\"node\" attr.name=\"kappa\" attr.type=\"double\" />"                   << std::endl;
  int k(9);
  for(int d(0); d<DIMENSION; ++d)
  {
    graphml_file << "  <key id=\"key" << k++ << "\" for=\"node\" attr.name=\"angle" << d << "\" attr.type=\"double\" />"  << std::endl;
  }
  if(OUTPUT_VERTICES_PROPERTIES)
  {
    graphml_file << "  <key id=\"key" << k++ << "\" for=\"node\" attr.name=\"expected_degree\" attr.type=\"double\" />"        << std::endl;
  }

  // Declares the graph.
  graphml_file << std::endl;
  graphml_file << "  <graph id=\"G\" edgedefault=\"undirected\" parse.nodeids=\"canonical\" parse.edgeids=\"canonical\" parse.order=\"nodesfirst\">" << std::endl;

  // Sets the graph properties.
  graphml_file << std::endl;
  graphml_file << "    <!-- graph properties -->" << std::endl;
  graphml_file << "    <data key=\"key0\">" << get_time()                << "</data>" << std::endl;
  graphml_file << "    <data key=\"key1\">" << HIDDEN_VARIABLES_FILENAME << "</data>" << std::endl;
  graphml_file << "    <data key=\"key2\">" << SEED                      << "</data>" << std::endl;
  graphml_file << "    <data key=\"key3\">" << DIMENSION                 << "</data>" << std::endl;
  graphml_file << "    <data key=\"key4\">" << BETA                      << "</data>" << std::endl;
  graphml_file << "    <data key=\"key5\">" << MU                        << "</data>" << std::endl;
  graphml_file << "    <data key=\"key6\">" << RADIUS                    << "</data>" << std::endl;

  // Declare the vertices
  graphml_file << std::endl;
  graphml_file << "    <!-- edges -->" << std::endl;
  for(int k, v(0); v<nb_vertices; ++v)
  {
    graphml_file << "    <node id=\"n" << v <<"\">" << std::endl;
    graphml_file << "      <data key=\"key7\">" << Num2Name[v] << "</data>"   << std::endl;
    graphml_file << "      <data key=\"key8\">" << kappa[v] << "</data>"   << std::endl;
    k = 9;
    for(int d(0); d<DIMENSION; ++d)
    {
      graphml_file << "      <data key=\"key" << k++ << "\">" << theta[v][d] << "</data>" << std::endl;

    }
    if(OUTPUT_VERTICES_PROPERTIES)
    {
      graphml_file << "      <data key=\"key" << k++ << "\">" << edegree[v] << "</data>" << std::endl;
    }
    graphml_file << "    </node>" << std::endl;
  }

  // Declare the edges
  int e = 0;
  graphml_file << std::endl;
  graphml_file << "    <!-- edges -->" << std::endl;
  for(auto edge : edgelist)
  {
    graphml_file << "    <edge id=\"e" << e <<"\" source=\"n" << edge.first <<"\" target=\"n" << edge.second << "\"/>" << std::endl;
    ++e;
  }

  // Footer.
  graphml_file << "  </graph>" << std::endl;
  graphml_file << "</graphml>" << std::endl;

  // Closes the stream.
  graphml_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::save_vertices_properties(std::vector<int>& rdegree, std::vector<double>& edegree, int width)
{
  // Sets the name of the file to write the hidden variables into.
  std::string vprop_filename = OUTPUT_ROOTNAME + "_vprop.dat";
  // Opens the stream and terminates if the operation did not succeed.
  std::fstream vprop_file(vprop_filename.c_str(), std::fstream::out);
  if( !vprop_file.is_open() )
  {
    std::cerr << "Could not open file: " << vprop_filename << "." << std::endl;
    std::terminate();
  }
  // Writes the header.
  vprop_file << "#";
  vprop_file << std::setw(width - 1) << "Vertex"         << " ";
  vprop_file << std::setw(width)     << "Kappa"          << " ";
  for(int d(0); d<DIMENSION; ++d)
  {
    vprop_file << std::setw(width)   << "Angle" + std::to_string(d + 1) << " ";
  }
  vprop_file << std::setw(width)     << "Act.Deg."       << " ";
  vprop_file << std::setw(width)     << "Exp.Deg."       << " ";
  vprop_file << std::endl;
  // Writes the hidden variables.
  for(int v(0); v<nb_vertices; ++v)
  {
    vprop_file << std::setw(width)   << Num2Name[v]      << " ";
    vprop_file << std::setw(width)   << kappa[v]         << " ";
    for(int d(0); d<DIMENSION; ++d)
    {
      vprop_file << std::setw(width) << theta[v][d]      << " ";
    }
    vprop_file << std::setw(width)   << rdegree[v]       << " ";
    vprop_file << std::setw(width)   << edegree[v]       << " ";
    vprop_file << std::endl;
  }
  // Closes the stream.
  vprop_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void modelSD_t::set_mu_to_default_value()
{
  // Makes sure the value of beta has been provided.
  if(BETA <= DIMENSION)
  {
    std::cerr << "ERROR: Default value for MU isn't valid if BETA <= DIMENSION" << std::endl;
    std::terminate();
  }
  // Computes the average value of kappa.
  double average_kappa = 0;
  for(int v(0); v<nb_vertices; ++v)
  {
    average_kappa += kappa[v];
  }
  average_kappa /= nb_vertices;
  // Sets MU to its default value.
  MU = BETA * std::sin(DIMENSION * PI / BETA) / (2.0 * PI * average_kappa);
  MU *= std::tgamma(static_cast<double>(DIMENSION) / 2.);
  MU /= std::pow(PI, static_cast<double>(DIMENSION) / 2.);
  // MU = BETA * std::sin(PI / BETA) / (2.0 * PI * average_kappa);
}





#endif // MODELSD_HPP_INCLUDED
