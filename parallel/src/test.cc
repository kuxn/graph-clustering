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
#include "test.h"
#include "graph.h"
#include "partition.h"
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;
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
    ifstream In("./test/par_test_8.dot");
    if (!In.is_open()) {
        std::cerr << "ERROR: Can't open the file" << endl;
        exit(-1);
    }
    Partition partition(g, 4, true);
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
    int num = 500, rank = 0;
    g.readDotFormat("./test/par_test_500.dot", num);
    if (world.rank() == rank) {
        g.printDotFormat();
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
        //g.printLaplacianMat();
    }
    env.~environment();
    return true;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testReadByColour
 *  Description:  Test reading the vertices with same colours and corresponding edges
 * =====================================================================================
 */

bool Tests::testReadByColour() {
    mpi::environment env;
    mpi::communicator world;
    Graph g;
    int num = 20, rank = 3;
    g.readDotFormatByColour("./test/par_test_20_4s.dot", num);
    if (world.rank() == rank) {
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
        g.printDotFormat();
    }
    env.~environment();
    return true;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testPartitionWithSubgraphs
 *  Description:  Cluster assignment: test the partition with partitioned graph
 * =====================================================================================
 */

bool Tests::testPartitionWithClusters() {
    mpi::environment env;
    mpi::communicator world;
    Graph g;
    int num = 1024, subgraphs = 4;
    bool gram_schmidt = true;
    g.readDotFormatByColour("./test/par_test_1024.dot", num);
    if (world.rank() == 0) {
        //g.printDotFormat();
        cout << "rank = " << world.rank() << ", size of graph = " << g.size() << ", num of edges = " << g.edgesNum() << endl;
    }
    Partition partition(g, subgraphs, gram_schmidt);
    env.~environment();
    return true;
}

// int main() {
//     Tests::testReadByColour();
//     //Tests::testPartition();
//     //Tests::testReadGraph();
//     //Tests::testPartitionWithClusters();
// 
//     return 0;
// }
