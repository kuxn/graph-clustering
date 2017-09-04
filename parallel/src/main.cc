/**
 * @file main.cc
 * @brief Main function
 * @author Ken Hu, xnchnhu@gmail.com
 */

#include <iostream>
#include <sstream>

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <boost/timer.hpp>

#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"
#ifdef VT_
#include "vt_user.h"
#endif

namespace mpi = boost::mpi;
namespace po = boost::program_options;
using namespace std;

int main(int argc, char* argv[]) {
#ifdef VT_
    VT_TRACER("MAIN");
#endif
    mpi::environment env;
    mpi::communicator world;

    int vertices, subgraphs;
    string filename;

    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", ":produce help message")
    ("output,o", ":output the partitioned graph into dot files")
    ("gram-schmidt,g", ":enable Gram Schmidt in Lanczos")
    ("read-by-colour,r", ":read dot format into different processes by colours")
    ("subgraphs,s", po::value<int>(), ":set number of subgraphs, has to be the power of 2")
    ("input-file,f", po::value<string>(), ":input file name")
    ("vertices,v", po::value<int>(), ":number of vertices of the input file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") && world.rank() == 0) {
        cout << desc << endl;
        return 0;
    }

    bool read_graph = vm.count("input-file") && vm.count("vertices"),
         read_by_colour = vm.count("read-by-colour"),
         sub_graphs = vm.count("subgraphs"),
         output = vm.count("output"),
         gram_schmidt = vm.count("gram-schmidt");

    Graph* g;

    boost::timer timer_io_input;
    if (read_graph) {
        g = new Graph;
        vertices = vm["vertices"].as<int>();
        filename = vm["input-file"].as<string>();
        if (world.rank() == 0) {
            cout << "Input file is \"" << filename  << "\""<< endl;
        }
        if (read_by_colour) {
            g->readDotFormatByColour(filename, vertices);
            if (world.rank() == 0) {
                cout << "cluster assignment" << endl;
            }
        } else {
            g->readDotFormat(filename, vertices);
            if (world.rank() == 0) {
                cout << "even assignment" << endl;
            }
        }
    } else {
        if (world.rank() == 0) {
            cout << "Please set file name and number of vertices" << endl;
            cout << desc << endl;
        }
        return 0;
    }
    if (sub_graphs) {
        subgraphs = vm["subgraphs"].as<int>();
        if (world.rank() == 0 && read_graph) {
            cout << "argument: subgraphs = " << subgraphs << "." << endl;
            cout << "number of processes = " << world.size() << "." << endl;
        }
    } else {
        subgraphs = 4;
        if (world.rank() == 0 && read_graph) {
            cout << "default argument: subgraphs = " << subgraphs << "." << endl;
            cout << "number of processes = " << world.size() << "." << endl;
        }
    }
    if (world.rank() == 0) {
        double t_input = timer_io_input.elapsed();
        cout << "input takes " << t_input << "s" << endl;
    }

    world.barrier();
    Partition partition(*g, subgraphs, gram_schmidt);
    world.barrier();

    boost::timer timer_io_output;
    if (world.rank() != 0) {
        filename = "./output/temp_";
        filename += to_string(vertices);
        filename += "v_";
        filename += to_string(world.size());
        filename += "wr_";
        filename += to_string(g->rank());
        filename += "r.dot";
        g->outputDotFormat(filename);
    }
    world.barrier();
    // Read all the vertices to rank 0
    if (world.rank() == 0) {
        for (int rank = 1; rank < world.size(); ++rank) {
            filename = "./output/temp_";
            filename += to_string(vertices);
            filename += "v_";
            filename += to_string(world.size());
            filename += "wr_";
            filename += to_string(rank);
            filename += "r.dot";
            g->readDotFormatWithColour(filename);
        }
        if (output) {
            filename = "./output/result_";
            filename += to_string(g->globalSize());
            filename += "v_";
            filename += to_string(subgraphs);
            filename += "s.dot";
            g->outputDotFormat(filename);
        }
        Analysis::outputTimes(world.size(), vertices, partition.times);
        double t_output = timer_io_output.elapsed();
        cout << "output takes " << t_output << "s" << endl;
        Analysis::cutEdgeVertexTable(*g, partition.ritz_values);
    }

    delete g;
    return 0;
}
