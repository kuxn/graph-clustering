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

	int num = 8;
	//cout << "num of vertices= " << num << endl;

    Graph g;
    ifstream In("test_8.dot");

	g.init(world.rank(), num, num/world.size());

    g.readDotFormat(In);

    if (world.rank() == rank) {
	    g.printDotFormat();
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
		g.printLaplacianMat();
    }

	

	return 0;
}
