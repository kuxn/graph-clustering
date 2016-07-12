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

#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"
#include "vt_user.h"

namespace mpi = boost::mpi;
namespace po = boost::program_options;
using namespace std;

int main(int argc, char* argv[]) {

    VT_TRACER("MAIN");
    mpi::environment env;
    mpi::communicator world;

    int vertices, subgraphs;
    string filename;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", ":produce help message")
        ("output,o", ":output the partitioned graph into dot files")
        ("gram-schmidt,g", ":Gram Schmidt in Lanczos")
        ("vertices,v", po::value<int>(), ":number of vertices of the input file")
        ("subgraphs,s", po::value<int>(), ":set number of subgraphs, has to be the power of 2")
        ("input-file,f", po::value<string>(), ":input file name")
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
        g->init(world.rank(), vertices, vertices/world.size());
        filename = vm["input-file"].as<string>();
        if (world.rank() == 0) {
            cout << "Input file is \"" << filename  << "\""<< endl;
        }
        g->readDotFormat(filename);
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
        }
    } else {
        subgraphs = 2;
        if (world.rank() == 0 && read_graph) {
            cout << "default argument: subgraphs = " << subgraphs << "." << endl;
        }
    }

    world.barrier();
    Partition partition(*g, subgraphs, gram_schmidt);
    world.barrier();

    if (output) {
        string filename("./output/parallel_");
        filename += to_string(g->localSize());
        filename += "v_";
        filename += to_string(g->rank());
        filename += "r.dot";
        g->outputDotFormat(filename);
    }

    world.barrier();
    // Read all the vertices to rank 0
    if (output && world.rank() == 0) {
        for (int rank = 0; rank < world.size(); rank++) {
            filename = "./output/parallel_";
            filename += to_string(vertices/world.size());
            filename += "v_";
            filename += to_string(rank);
            filename += "r.dot";
            g->readDotFormatWithColour(filename);
        }
        g->printDotFormat();
        partition.printLapEigenvalues();
        partition.printLapEigenMat();
        Analysis::cutEdgeVertexTable(*g);
    }
    //fstream Output;
    //for (int proc = 0; proc < world.size(); proc++) {
    //	if (world.rank() == proc) {
    //		//Output.open("output.dot", std::ofstream::ate);
    //		Output.open("output.dot", std::fstream::in | std::fstream::out | std::fstream:: ate);
    //		for (int vertex = 0; vertex < g.localSize(); vertex++) {
    //			//Output << g.globalIndex(vertex) << "[Colour=" << g.getColour(g.globalIndex(vertex)) << "];" << endl;	
    //			cout << g.globalIndex(vertex) << "[Colour=" << g.getColour(g.globalIndex(vertex)) << "];" << endl;	
    //		}
    //		Output.close();
    //	}
    //}
    //world.barrier();
    //for (int proc = 0; proc < world.size(); proc++) {
    //	if (world.rank() == proc) {
    //		//Output.open("output.dot", std::ofstream::ate);
    //		Output.open("output.dot", std::fstream::in | std::fstream::out | std::fstream:: ate);
    //		for (int vertex = 0; vertex < g.localSize(); vertex++) {
    //			auto it = g.find(g.globalIndex(vertex));
    //			for (auto neighbour:it->second)
    //				//Output << g.globalIndex(vertex) << "--" << neighbour << ";" << endl;
    //				cout << g.globalIndex(vertex) << "--" << neighbour << ";" << endl;
    //		}
    //		Output.close();
    //	}
    //}

    //env.~environment();
    //cout << "finalized = " << env.finalized() << endl;

    delete g;
    return 0;
}
