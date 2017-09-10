/**
 * @file main.cc
 * @brief Main function
 * @author Ken Hu, xnchnhu@gmail.com
 */

#include <boost/program_options.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include <string>

#include "analysis.h"
#include "graph.h"
#include "partition.h"

#ifdef VT_
#include "vt_user.h"
#endif

using namespace std;
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
#ifdef VT_
    VT_TRACER("MAIN");
#endif
    int vertices, colours;
    string filename;

    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", ":produce help message")
    ("output,o", ":output the partitioned graph into a dot file")
    ("gram-schmidt,g", ":enable Gram Schmidt in Lanczos, default: false")
    ("benchmarks,b", ":run benchmarks, default: false")
    ("vertices,v", po::value<int>(), ":set number of vertices, default: 20")
    ("colours,c", po::value<int>(), ":set number of colours, has to be the power of 2, default: 2")
    ("input-file,f", po::value<string>(), ":input file name")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") || vm.size() == 0) {
        cout << desc << endl;
        return 1;
    }

    bool random_graph = vm.count("vertices"),
         read_graph = vm.count("input-file"),
         sub_graphs = vm.count("colours"),
         output = vm.count("output"),
         gram_schmidt = vm.count("gram-schmidt"),
         benchmarks = vm.count("benchmarks");

    Graph* g;

    if (random_graph && read_graph) {
        cout << "WARNING: you can either set vertices or read graphs" << endl;
        cout << desc << endl;
        return 1;
    }
    if (random_graph) {
        vertices = vm["vertices"].as<int>();
        cout << "argument: vertices =  " << vertices << "." << endl;
        g = new Graph(vertices);
    } else if (read_graph) {
        g = new Graph;
        filename = vm["input-file"].as<string>();
        cout << "Input file is \"" << filename << "\"" << endl;
        g->readDotFormat(filename);
    } else {
        cout << desc << endl;
        vertices = 20;
        cout << "default argument: vertices = " << vertices << "." << endl;
        g = new Graph(vertices);
    }
    if (sub_graphs) {
        colours = vm["colours"].as<int>();
        cout << "argument: colours = " << colours << "." << endl;
    } else {
        colours = 2;
        cout << "default argument: colours = " << colours << "." << endl;
    }

    if (benchmarks) {
        Analysis::benchmarks(gram_schmidt);
    } else {
        Partition partition(*g, colours, gram_schmidt);
        if (output) {
            string filename("./output/serial_");
            filename += to_string(g->size());
            filename += "v_";
            filename += to_string(g->subgraphsNum());
            filename += "s.dot";
            g->outputDotFormat(filename);
        }
        // g->outputDotFormat("plot_6.dot");
        // g->printLaplacianMat();
        // partition.printLapEigenvalues();
        // partition.outputLapEigenvalues();
        // partition.printLapEigenMat();
        Analysis::outputTimes(g->size(), partition.times);
        Analysis::cutEdgeVertexTable(*g, partition.ritzValues);
    }

    delete g;
    return 0;
}
