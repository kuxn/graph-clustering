/**
 * @file ParallelTest.cc
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
    virtual void SetUp()
    {
        char* path = std::getenv("TEST_FILES");
        if (path != nullptr) {
            filePath = path;
        } else {
            filePath = "../../tests/dotfiles";
            // std::cout << "Please export TEST_FILES to the path of test files"
            //<< std::endl;
        }
    }
    Graph g;
    std::string filePath;
    mpi::communicator world;
};

/**
 * @brief Partitioning is based on the correctness of the calculation of
 *        eigenvalues and corresponding eigenvectors of the Laplacian matrix
 */
TEST_F(ParallelTest, ReadDotFilReadDotFilee)
{
    int num = 500;
    int rank = 0;
    g.readDotFormat(filePath + "/par_test_500.dot", num);
    if (world.rank() == rank) {
        // g.printDotFormat();
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
        // g.printLaplacianMat();
    }
    ASSERT_TRUE(true);
}

TEST_F(ParallelTest, testPartition)
{
    int num = 8;
    int rank = 0;
    g.readDotFormat(filePath + "/par_test_8.dot", num);
    Partition partition(g, 4, true);
    if (world.rank() == rank) {
        partition.printLapEigenvalues();
        partition.printLapEigenMat();
        // g.printDotFormat();
    }
}

/**
 * @brief Test the read function, can not verify the correctness in unit
 *        testing, but can test the performance
 */
TEST_F(ParallelTest, testReadGraph)
{
    int num = 500;
    int rank = 0;
    g.readDotFormat(filePath + "/par_test_500.dot", num);
    if (world.rank() == rank) {
        // g.printDotFormat();
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
        // g.printLaplacianMat();
    }
}

/**
 * @brief Test reading the vertices with same colours and corresponding edges
 */
TEST_F(ParallelTest, testReadGraphByColour)
{
    int num = 20;
    int rank = 3;
    g.readDotFormatByColour(filePath + "/par_test_20_4s.dot", num);
    if (world.rank() == rank) {
        cout << "rank = " << world.rank() << endl;
        cout << "size of graph = " << g.size() << endl;
        cout << "num of edges = " << g.edgesNum() << endl;
        g.printDotFormat();
    }
}

/**
 * @brief Cluster assignment: test the partition with partitioned graph
 */
TEST_F(ParallelTest, testPartitionWithClusters)
{
    int num = 1024;
    int subgraphs = 4;
    bool gram_schmidt = true;
    g.readDotFormatByColour(filePath + "/par_test_1024_4s.dot", num);
    if (world.rank() == 0) {
        // g.printDotFormat();
        cout << "rank = " << world.rank() << ", size of graph = " << g.size()
             << ", num of edges = " << g.edgesNum() << endl;
    }
    Partition partition(g, subgraphs, gram_schmidt);
}

int main(int argc, char** argv)
{
    int result = 0;
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env;
    mpi::communicator world;
    // Gets hold of the event listener list.
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    // Only leave the listener for rank 0, such that the messages don't mingle
    if (world.rank() != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }

    result = RUN_ALL_TESTS();
    env.~environment();
    return result;
}
