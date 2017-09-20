/**
 * @file SerialTest.cc
 * @brief Functions for unit tests
 * @author Ken Hu, xnchnhu@gmail.com
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include "analysis.h"
#include "graph.h"
#include "gtest/gtest.h"
#include "lanczos.h"
#include "partition.h"
#include "tqli.h"

using namespace std;

class SerialTest : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        char* path = std::getenv("TEST_FILES");
        if (path != nullptr) {
            filePath = path;
        } else {
            std::cout << "Please export TEST_FILES to the path of test files"
                      << std::endl;
        }
    }
    Graph g;
    std::string filePath;
};

/**
 * @brief Test the read function, can not verify the correctness in unit
 *        testing, but can test the performance
 */
TEST_F(SerialTest, testReadGraph)
{
    g.readDotFormat(filePath + "/test_read_20.dot");

    EXPECT_EQ(20, g.size());
    EXPECT_EQ(36, g.edgesNum());
    EXPECT_EQ(1, g.subgraphsNum());
}

TEST_F(SerialTest, testReadGraphWithColour)
{
    g.readDotFormatWithColour(filePath + "/test_read_20.dot");

    EXPECT_EQ(g.subgraphsNum(), 4);
}

/**
 * @brief Test tqli function for the correctness of calculating eigenvalues,
 *        correct eigenvalues come from Matlab.
 */
TEST_F(SerialTest, testTqli)
{
    int size = 5;
    vector<vector<double>> eigenvecs;
    vector<double> diagonal, subdiagonal;
    diagonal = {0.569893, 3.81259, 3.02478, 3.39064, 3.2021};
    subdiagonal = {1.45159, 0.550477, 1.06987, 1.25114, 0.0};

    tqli(diagonal, subdiagonal, eigenvecs);
    vector<double> result = {0.0, 1.58578, 3.00000, 4.41421, 5.00000};

    for (int i = 0; i < size; i++) {
        EXPECT_LE(diagonal[i] - result[i], 1e-5);
    }
}

/**
 * @brief The output alpha/beta would vary based on different initial vector,
 *        the eigenvalues should always be same.
 */
TEST_F(SerialTest, testLanczos)
{
    Graph g;
    g.readDotFormat(filePath + "/par_test_8.dot");
    int size = g.size();

    // Calculate the diagonal and subdiagonal vectors
    Lanczos<vector<double>, double> lanczos(g, size, true);
    vector<double> alpha = lanczos.alpha;
    vector<double> beta = lanczos.beta;

    beta.push_back(0);

    // Create the identity matrix used as input for TQLI
    vector<vector<double>> eigenvecs;
    tqli(alpha, beta, eigenvecs);
    vector<double> eigenvalues = {0,       1.20972, 1.505, 2,
                                  2.86246, 4.32623, 5,     7.09659};

    sort(alpha.begin(), alpha.end());
    for (int i = 0; i < size; i++) {
        EXPECT_LT(abs(alpha[i] - eigenvalues[i]), 1e-5);
    }
}

/**
 * @brief Partitioning is based on the correctness of the calculation of
 *        eigenvalues and corresponding eigenvectors of the Laplacian matrix
 */
TEST_F(SerialTest, testPartition)
{
    Graph g;
    g.readDotFormat(filePath + "/test_partition_10.dot");

    Partition partition(g, 2, true);
    double cut_edge_percent = Analysis::cutEdgePercent(g);

    EXPECT_EQ(g.subgraphsNum(), 2);
    EXPECT_LT(std::abs(cut_edge_percent - 0.142857), 1e-5);
}

/**
 * @brief Test the function manuallyPartition in analysis.cc
 */
TEST_F(SerialTest, testRandomPartition)
{
    g.readDotFormat(filePath + "/test_1000.dot");
    int colours = 4;
    Analysis::randomPartition(g, colours);
    double cut_edge_percent = Analysis::cutEdgePercent(g);

    EXPECT_EQ(g.subgraphsNum(), 4);
    EXPECT_LT(std::abs(cut_edge_percent - 0.740337), 1e-5);
}

TEST_F(SerialTest, testEvenPartition)
{
    g.readDotFormat(filePath + "/test_1000.dot");
    int colours = 4;
    Analysis::evenPartition(g, colours);
    double cut_edge_percent = Analysis::cutEdgePercent(g);
    // std::vector<double> ritz_values(1, 0.0);
    // Analysis::cutEdgeVertexTable(g, ritz_values);

    EXPECT_EQ(g.subgraphsNum(), 4);
    EXPECT_LT(std::abs(cut_edge_percent - 0.751239), 1e-5);
}

/**
 * @brief Test the connection of the graph after partitioning
 */
TEST_F(SerialTest, testCutEdgeVertexTable)
{
    g.readDotFormatWithColour(filePath + "/test_read_20.dot");
    vector<double> ritz_value = {0.868758, 1.1268};
    // Analysis::cutEdgeVertexTable(g, ritz_value);

    /*-----------------------------------------------------------------------------
     * Basic info of the graph
     *-----------------------------------------------------------------------------
    Vertices:   20
    Edges:      36
    Subgraphs:  4
    Used Ritz values: 0.868758, 1.1268
    Cut Edge Percent: 47.2222%
     *-----------------------------------------------------------------------------
     * Number of nodes in each subgraph
     *-----------------------------------------------------------------------------
    Colour:     0      1       2       3
    Vertices:   6      5       3       6
     *-----------------------------------------------------------------------------
     * Edges table after partitioning
     * Each element represents number of edges (inside)/between subgraphs
     *-----------------------------------------------------------------------------
        0      1      2      3
     0 (4)     4      5      3
     1 4      (6)     0      3
     2 5       0     (2)     2
     3 3       3      2     (7)
     *-----------------------------------------------------------------------------*/
}

TEST_F(SerialTest, testReothogonalisation)
{
    Graph g(100);
    Partition partition1(g, 4, false);
    // cout << "WITHOUT reorthogonalisation: " << endl;
    // Analysis::cutEdgeVertexTable(g, partition1.ritzValues);
    // cout << "eigenvalues:";
    // partition1.printLapEigenvalues();
    // partition1.printLapEigenMat();

    Partition partition2(g, 4, true);
    // cout << "WITH reorthogonalisation: " << endl;
    // Analysis::cutEdgeVertexTable(g, partition2.ritzValues);
    // cout << "eigenvalues:";
    // partition2.printLapEigenvalues();
    // partition2.printLapEigenMat();
}
