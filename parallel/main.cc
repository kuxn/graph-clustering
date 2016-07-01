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
#include <boost/optional/optional_io.hpp>

#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"

namespace mpi = boost::mpi;
using namespace std;

int main(int argc, char* argv[]) {

    mpi::environment env;
    mpi::communicator world;

    int rank;
	istringstream ss(argv[1]);
	if (!(ss >> rank))
		cerr << "Invalid number " << argv[1] << endl;

	int num = 200;
	//cout << "num of vertices= " << num << endl;

    Graph g;
    //ifstream In("par_test_8.dot");
    //ifstream In("par_test_500.dot");
    ifstream In("par_test_200.dot");

	g.init(world.rank(), num, num/world.size());

    g.readDotFormat(In);

    //if (world.rank() == rank) {
	//    g.printDotFormat();
    //    cout << "rank = " << world.rank() << endl;
    //    cout << "size of graph = " << g.size() << endl;
    //    cout << "num of edges = " << g.edgesNum() << endl;
	//	g.printLaplacianMat();
    //}

	world.barrier();

    Partition partition(world, g, 8, true);
	if (world.rank() == 0) {
		partition.printLapEigenvalues();
		partition.printLapEigenMat();
		//g.printDotFormat();
	}
	
	world.barrier();
	// Print coloured graph from each process to a single file

	string filename("output_rank_");
	filename += to_string(world.rank());
	filename += ".dot";

	g.outputDotFormat(filename);

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

	env.~environment();
	//cout << "finalized = " << env.finalized() << endl;

	return 0;
}

