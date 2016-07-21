/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  Main function
 *        Created:  06/09/2016 18:12:52
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/timer.hpp>

#include "graph.h"
#include "partition.h"
#include "analysis.h"
#include "vt_user.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char* argv[]) {

    VT_TRACER("MAIN");
    int vertices, subgraphs;
    string filename;

    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", ":produce help message")
    ("output,o", ":output the partitioned graph into a dot file")
    ("gram-schmidt,g", ":Enable Gram Schmidt in Lanczos, default: false")
    ("vertices,v", po::value<int>(), ":set number of vertices, default: 20")
    ("subgraphs,s", po::value<int>(), ":set number of subgraphs, has to be the power of 2, default: 2")
    ("input-file,f", po::value<string>(), ":input file name")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    bool random_graph = vm.count("vertices"),
         read_graph = vm.count("input-file"),
         sub_graphs = vm.count("subgraphs"),
         output = vm.count("output"),
         gram_schmidt = vm.count("gram-schmidt");

    Graph* g;

    if (random_graph && read_graph) {
        cout << "WARNING: you can either set vertices or read graphs" << endl;
        cout << desc << "\n";
        return 1;
    }
    if (random_graph) {
        vertices = vm["vertices"].as<int>();
        cout << "argument: vertices =  " << vertices << "." << endl;
        g = new Graph(vertices);
    } else if (read_graph) {
        g = new Graph;
        filename = vm["input-file"].as<string>();
        cout << "Input file is \"" << filename  << "\""<< endl;
        g->readDotFormat(filename);
    } else {
        cout << desc << "\n";
        vertices = 20;
        cout << "default argument: vertices = " << vertices << "." << endl;
        g = new Graph(vertices);
    }
    if (sub_graphs) {
        subgraphs = vm["subgraphs"].as<int>();
        cout << "argument: subgraphs = " << subgraphs << "." << endl;
    } else {
        subgraphs = 2;
        cout << "default argument: subgraphs = " << subgraphs << "." << endl;
    }

	boost::timer timer_partition;
    Partition partition(*g, subgraphs, gram_schmidt);
	cout << "Partition takes " << timer_partition.elapsed() << "s" << endl;

    if (output) {
        string filename("./output/serial_");
        filename += to_string(g->size());
        filename += "v_";
        filename += to_string(g->subgraphsNum());
        filename += "s.dot";
        g->outputDotFormat(filename);
    } else {
        //g->printDotFormat();
        //g->printLaplacianMat();
        partition.printLapEigenvalues();
        //partition.printLapEigenMat();
        Analysis::cutEdgeVertexTable(*g);
    }

    delete g;
    return 0;
}
