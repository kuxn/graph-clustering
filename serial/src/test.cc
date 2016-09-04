/*
 * =====================================================================================
 *
 *       Filename:  test.cc
 *
 *    Description:  Functions for unit tests
 *        Created:  06/09/2016 18:12:52
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <algorithm>
#include "test.h"
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"

using namespace std;

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testReadGraph
 *  Description:  Test the read function, can not verify the correctness in unit testing, but can test the performance
 * =====================================================================================
 */

bool Tests::testReadGraph() {
    Graph g;
    g.readDotFormat("./test/test_read_20.dot");
    if (g.size() != 20 || g.edgesNum() != 36 || g.subgraphsNum() != 1) {
        return false;
    }
    return true;
}

bool Tests::testReadGraphWithColour() {
    Graph g;
    g.readDotFormatWithColour("./test/test_read_20.dot");
    if (g.subgraphsNum() != 4) {
        return false;
    }
    return true;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testTqli
 *  Description:  Test tqli function for the correctness of calculating eigenvalues, correct eigenvalues come from Matlab.
 * =====================================================================================
 */

bool Tests::testTqli() {
    int size = 5;
    vector<vector<double>> eigenvecs;
    vector<double> diagonal, subdiagonal;
    diagonal = {0.569893, 3.81259, 3.02478, 3.39064, 3.2021};
    subdiagonal = {1.45159, 0.550477, 1.06987, 1.25114, 0.0};

    tqli(diagonal, subdiagonal, eigenvecs);
    vector<double> result = {0.0, 1.58578, 3.00000, 4.41421, 5.00000};

    for (int i = 0; i < size; i++) {
        if (std::abs(diagonal[i] - result[i]) > 1e-5)
            return false;
    }
    return true;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testLanczos
 *  Description:  The output alpha/beta would vary based on different initial vector, the eigenvalues should always be same.
 * =====================================================================================
 */

bool Tests::testLanczos() {
    Graph g;
    g.readDotFormat("./test/test_lanczos_8.dot");
    int size = g.size();

    // Calculate the diagonal and subdiagonal vectors
    Lanczos<vector<double>, double> lanczos(g, size, true);
    vector<double> alpha = lanczos.alpha;
    vector<double> beta = lanczos.beta;

    beta.push_back(0);

    // Create the identity matrix used as input for TQLI
    vector<vector<double>> eigenvecs;
    tqli(alpha, beta, eigenvecs);
    vector<double> eigenvalues = {0, 1.20972, 1.505, 2, 2.86246, 4.32623, 5, 7.09659};

    sort(alpha.begin(), alpha.end());
    for (int i = 0; i < size; i++) {
        if (abs(alpha[i] - eigenvalues[i]) > 1e-5)
            return false;
    }
    return true;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testPartition
 *  Description:  Partitioning is based on the correctness of the calculation of eigenvalues and corresponding eigenvectors of the Laplacian matrix
 * =====================================================================================
 */

bool Tests::testPartition() {
    Graph g;
    g.readDotFormat("./test/test_partition_10.dot");

    Partition partition(g, 2, true);
    double cut_edge_percent = Analysis::cutEdgePercent(g);

    if (g.subgraphsNum() != 2 || std::abs(cut_edge_percent - 0.142857) > 1e-5) {
        return false;
    }
    return true;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testManuallyPartition
 *  Description:  Test the function manuallyPartition in analysis.cc
 * =====================================================================================
 */

bool Tests::testRandomPartition() {
    Graph g;
    g.readDotFormat("./test/test_1000.dot");
    int colours = 4;
    Analysis::randomPartition(g, colours);
    double cut_edge_percent = Analysis::cutEdgePercent(g);
    cout << "cut_edge_percent = " << cut_edge_percent * 100 << "%"<< endl;

    if (g.subgraphsNum() != 4 || std::abs(cut_edge_percent - 0.501239) < 1e-5) {
        return false;
    }
    return true;
}

bool Tests::testEvenPartition() {
    Graph g;
    g.readDotFormat("./test/test_1000.dot");
    int colours = 4;
    Analysis::evenPartition(g, colours);
    double cut_edge_percent = Analysis::cutEdgePercent(g);
    std::vector<double> ritz_values(1, 0.0);
    Analysis::cutEdgeVertexTable(g, ritz_values);
    cout << "cut_edge_percent = " << cut_edge_percent * 100 << "%"<< endl;

    if (g.subgraphsNum() != 4 || std::abs(cut_edge_percent - 0.501239) < 1e-5) {
        return false;
    }
    return true;
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  testCutEdgeTable
 *  Description:  Test the connection of the graph after partitioning
 * =====================================================================================
 */

bool Tests::testCutEdgeVertexTable() {
    Graph g;

    g.readDotFormatWithColour("./test/test_read_20.dot");
    vector<double> ritz_value = {0.868758, 1.1268};
    Analysis::cutEdgeVertexTable(g, ritz_value);

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
    Colour: 	0	1	2	3
    Vertices: 	6	5	3	6
     *-----------------------------------------------------------------------------
     * Edges table after partitioning
     * Each element represents number of edges (inside)/between subgraphs
     *-----------------------------------------------------------------------------
     	0	1	2	3
     0	(4)	4	5	3
     1	4	(6)	0	3
     2	5	0	(2)	2
     3	3	3	2	(7)
     *-----------------------------------------------------------------------------*/
    return true;
}

bool Tests::testReothogonalisation() {
    Graph g(100);
    Partition partition1(g, 4, false);
    cout << "WITHOUT reorthogonalisation: " << endl;
    Analysis::cutEdgeVertexTable(g, partition1.ritz_values);
    cout << "eigenvalues:";
    //partition1.printLapEigenvalues();
    //partition1.printLapEigenMat();

    Partition partition2(g, 4, true);
    cout << "WITH reorthogonalisation: " << endl;
    Analysis::cutEdgeVertexTable(g, partition2.ritz_values);
    //cout << "eigenvalues:";
    //partition2.printLapEigenvalues();
    //partition2.printLapEigenMat();

    return true;
}

//int main() {
    //Tests::testLanczos();
    //Tests::testCutEdgeVertexTable();
    //Tests::testRandomPartition();
    //Tests::testEvenPartition();
    //Tests::testReothogonalisation();

    //return 0;
//}
