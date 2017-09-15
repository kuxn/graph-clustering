/**
 * @file test.cc
 * @brief Functions for unit tests
 * @author Ken Hu, xnchnhu@gmail.com
 */

#include <boost/mpi.hpp>
#include <fstream>
#include <iostream>
#include <utility>
#include "graph.h"
#include "gtest/gtest.h"
#include "partition.h"

namespace mpi = boost::mpi;
using namespace std;

class ParallelTest : public ::testing::Test
{
protected:
    virtual void SetUp() {}
};

/**
 * @brief Partitioning is based on the correctness of the calculation of
 * eigenvalues and corresponding eigenvectors of the Laplacian matrix
 * @param FILL-ME-IN
 * @return FILL-ME-IN
 */
TEST(ParallelTest, ReadDotFile) {
    mpi::communicator world;
    Graph g;
    int num = 500, rank = 0;
    g.readDotFormat("../tests/dotfiles/par_test_500.dot", num);
    if (world.rank() == rank) {
        //g.printDotFormat();
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
        // g.printLaplacianMat();
    }
    ASSERT_TRUE(true);
}

int main(int argc, char** argv) {
    int result = 0;
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env;
    mpi::communicator world;
    // Gets hold of the event listener list.
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    // Only leave the listener for rank 0, such that the messages dorn't mingle
    if (world.rank() != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }

    return result;
}

//bool Tests::testPartition()
//{
//    mpi::environment env;
//    mpi::communicator world;
//    Graph g;
//    ifstream In("./test/par_test_8.dot");
//    if (!In.is_open()) {
//        std::cerr << "ERROR: Can't open the file" << endl;
//        exit(-1);
//    }
//    Partition partition(g, 4, true);
//    if (world.rank() == 0) {
//        partition.printLapEigenvalues();
//        partition.printLapEigenMat();
//        g.printDotFormat();
//    }
//    env.~environment();
//    return true;
//}
//
///**
// * @brief Test the read function, can not verify the correctness in unit
// * testing, but can test the performance
// * @param FILL-ME-IN
// * @return FILL-ME-IN
// */
//
//bool Tests::testReadGraph()
//{
//    mpi::environment env;
//    mpi::communicator world;
//    Graph g;
//    int num = 500, rank = 0;
//    g.readDotFormat("./test/par_test_500.dot", num);
//    if (world.rank() == rank) {
//        g.printDotFormat();
//        cout << "rank = " << world.rank() << endl;
//        cout << "size of graph = " << g.size() << endl;
//        cout << "num of edges = " << g.edgesNum() << endl;
//        // g.printLaplacianMat();
//    }
//    env.~environment();
//    return true;
//}
//
///**
// * @brief Test reading the vertices with same colours and corresponding edges
// * @param FILL-ME-IN
// * @return FILL-ME-IN
// */
//
//bool Tests::testReadByColour()
//{
//    mpi::environment env;
//    mpi::communicator world;
//    Graph g;
//    int num = 20, rank = 3;
//    g.readDotFormatByColour("./test/par_test_20_4s.dot", num);
//    if (world.rank() == rank) {
//        cout << "rank = " << world.rank() << endl;
//        cout << "size of graph = " << g.size() << endl;
//        cout << "num of edges = " << g.edgesNum() << endl;
//        g.printDotFormat();
//    }
//    env.~environment();
//    return true;
//}
//
///**
// * @brief Cluster assignment: test the partition with partitioned graph
// * @param FILL-ME-IN
// * @return FILL-ME-IN
// */
//
//bool Tests::testPartitionWithClusters()
//{
//    mpi::environment env;
//    mpi::communicator world;
//    Graph g;
//    int num = 1024, subgraphs = 4;
//    bool gram_schmidt = true;
//    g.readDotFormatByColour("./test/par_test_1024.dot", num);
//    if (world.rank() == 0) {
//        // g.printDotFormat();
//        cout << "rank = " << world.rank() << ", size of graph = " << g.size()
//             << ", num of edges = " << g.edgesNum() << endl;
//    }
//    Partition partition(g, subgraphs, gram_schmidt);
//    env.~environment();
//    return true;
//}

