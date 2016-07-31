/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  Main function
 *        Created:  06/29/2016 18:12:52
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <sstream>

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <boost/mpi/timer.hpp>

#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"
//#include "vt_user.h"

namespace mpi = boost::mpi;
namespace po = boost::program_options;
using namespace std;

int main(int argc, char* argv[]) {

    //VT_TRACER("MAIN");
    mpi::environment env;
    mpi::communicator world;

    int vertices, subgraphs;
    string filename;

    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", ":produce help message")
    ("output,o", ":output the partitioned graph into dot files")
    ("gram-schmidt,g", ":Gram Schmidt in Lanczos")
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
         sub_graphs = vm.count("subgraphs"),
         output = vm.count("output"),
         gram_schmidt = vm.count("gram-schmidt");

    Graph* g;

    if (read_graph) {
        g = new Graph;
        vertices = vm["vertices"].as<int>();
        filename = vm["input-file"].as<string>();
        if (world.rank() == 0) {
            cout << "Input file is \"" << filename  << "\""<< endl;
        }
        g->readDotFormat(filename, vertices);
        //g->readDotFormatByColour(filename, vertices);
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

    world.barrier();
    Partition partition(*g, subgraphs, gram_schmidt);
    std::vector<int> local_size_lookup;
    all_gather(world, g->localSize(), local_size_lookup);
    world.barrier();

    if (output && world.rank() != 0) {
        filename = "./output/temp_";
        filename += to_string(vertices);
        filename += "v_";
        filename += to_string(g->localSize());
        filename += "lv_";
        filename += to_string(g->rank());
        filename += "r.dot";
        g->outputDotFormat(filename);
    }
    world.barrier();
    // Read all the vertices to rank 0
    if (world.rank() == 0) {
        if (output) {
            for (int rank = 1; rank < world.size(); rank++) {
                filename = "./output/temp_";
                filename += to_string(vertices);
                filename += "v_";
                filename += to_string(local_size_lookup[rank]);
                filename += "lv_";
                filename += to_string(rank);
                filename += "r.dot";
                g->readDotFormatWithColour(filename);
            }
            filename = "./output/result_";
            filename += to_string(g->globalSize());
            filename += "v_";
            filename += to_string(subgraphs);
            filename += "s.dot";
            g->outputResult(filename);
        }
        //g->printDotFormat();
        partition.printLapEigenvalues();
        //partition.printLapEigenMat();
        Analysis::outputTimes(world.size(), vertices, partition.times);
        Analysis::cutEdgeVertexTable(*g, partition.ritz_values);
    }
	//partition.printLapEigenMat();

    env.~environment();

    delete g;
    return 0;
}
