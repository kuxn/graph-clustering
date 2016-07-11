/*
 * =====================================================================================
 *
 *       Filename:  test.cc
 *
 *    Description:  Functions for unit tests
 *        Created:  06/29/2016 18:12:52
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include "test.h"
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"

using namespace std;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testPartition
 *  Description:  Partitioning is based on the correctness of the calculation of eigenvalues and corresponding eigenvectors of the Laplacian matrix
 * =====================================================================================
 */

bool Tests::testPartition() {

    mpi::environment env;
    mpi::communicator world;

    Graph g;
    ifstream In("par_test_8.dot");
    if (!In.is_open()) {
        std::cerr << "ERROR: Can't open the file" << endl;
        exit(-1);
    }
    Partition partition(world, g, 4, true);
    if (world.rank() == 0) {
        partition.printLapEigenvalues();
        partition.printLapEigenMat();
        g.printDotFormat();
    }

    env.~environment();

    return true;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testReadGraph
 *  Description:  Test the read function, can not verify the correctness in unit testing, but can test the performance
 * =====================================================================================
 */

bool Tests::testReadGraph() {

    mpi::environment env;
    mpi::communicator world;

    Graph g;
    //ifstream In("par_test_8.dot");
    ifstream In("par_test_500.dot");
    if (!In.is_open()) {
        std::cerr << "ERROR: Can't open the file" << endl;
        exit(-1);
    }

    g.init(world.rank(), num, num/world.size());
    g.readDotFormat(In);

    int rank = 0;
    if (world.rank() == rank) {
        g.printDotFormat();
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
        g.printLaplacianMat();
    }

    env.~environment();
    return true;
}

int main() {

    Tests::testReadGraph();
    //Tests::testPartition();

    return 0;
}
