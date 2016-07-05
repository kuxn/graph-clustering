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
		("vertices,v", po::value<int>(), ":set number of vertices")
		("subgraphs,s", po::value<int>(), ":set number of subgraphs, has to be the power of 2")
		("input-file,f", po::value<string>(), ":input file name")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}
	if (vm.count("subgraphs")) {
		subgraphs = vm["subgraphs"].as<int>();
		cout << "argument: subgraphs = " << subgraphs << "." << endl;	
	} else {
		subgraphs = 2;
		cout << "default argument: subgraphs = " << subgraphs << "." << endl;	
	}

	bool random_graph = vm.count("vertices"), 
	     read_graph = vm.count("input-file"), 
	     output = vm.count("output");

	Graph* g;

	if (random_graph) {
		vertices = vm["vertices"].as<int>();
		cout << "argument: vertices =  " << vertices << "." << endl;	
		g = new Graph(vertices);
	} else if (read_graph) {
		g = new Graph;
		filename = vm["input-file"].as<string>();
		cout << "Input file is \"" << filename  << "\""<< endl;
		ifstream In(filename);
		g->readDotFormat(In);
	} else {
		vertices = 20;
		cout << "default argument: vertices = " << vertices << "." << endl;	
		g = new Graph(vertices);
	}
	
	Partition partition(*g, subgraphs, true);

	if (output) {
		string filename("serial_");
		filename += to_string(g->size());
		filename += "v_";
		filename += to_string(g->subgraphsNum());
		filename += "s.dot";
		g->outputDotFormat(filename);
	} else {
		//g->printDotFormat();
		//g->printLaplacianMat();
		//partition.printLapEigenvalues();
		//partition.printLapEigenMat();
		Analysis::cutEdgeVertexTable(*g);
	}

	delete g;
	return 0;
}

