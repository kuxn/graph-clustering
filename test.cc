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
#include "test.h"
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"

namespace Tests {

using namespace std;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testTqli
 *  Description:  Test tqli function for the correctness of calculating eigenvalues, correct eigenvalues come from Matlab.
 * =====================================================================================
 */

bool testTqli() {
	int size = 5;
	unordered_map<int, vector<double>> eigenvecs;
	vector<double> vinitial(size, 0);
	for(int i = 0; i < size; i++) {
		eigenvecs.insert({i, vinitial});
	}

	for(int i = 0; i < size; i++) {
		eigenvecs[i][i] = 1;
	}

	vector<double> diagonal, subdiagonal;
	diagonal = {0.569893, 3.81259, 3.02478, 3.39064, 3.2021};
	subdiagonal = {1.45159, 0.550477, 1.06987, 1.25114, 0.0};

	tqli(diagonal, subdiagonal, size, eigenvecs);
	vector<double> result = {0.0, 1.58578, 3.00000, 4.41421, 5.00000};
	
	for (int i = 0; i < size; i++) {
		if (abs(diagonal[i] - result[i]) > 1e-5)
			return false;
	}
	
	return true;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testLanczos
 *  Description:  The eigenvalues of original matrix should be same as the ones after tridiagonalising and calculating by TQLI.
 * =====================================================================================
 */

bool testLanczos() {
	// Initialise the graph
	Graph g;
	g.addEdge(0,1);
	g.addEdge(0,4);
	g.addEdge(1,2);
	g.addEdge(2,3);
	g.addEdge(2,4);
	g.addEdge(3,4);

	// The Laplacian matrix of the undirected graph above:
	/*-----------------------------------------------------------------------------
			0	1	2	3	4
		0	2	-1	0	0	-1
		1	-1	2	-1	0	0
		2	0	-1	3	-1	-1
		3	0	0	-1	2	-1
		4	-1	0	-1	-1	3
	 *-----------------------------------------------------------------------------*/

	int size = g.size();

	// Calculate the diagonal and subdiagonal vectors
	Lanczos<vector<double>, double> lanczos(g, true);
	vector<double> alpha = lanczos.alpha;
	vector<double> beta = lanczos.beta;

	beta.push_back(0);

	// Initialise the input matrix for storing eigenvectors
	unordered_map<int, vector<double>> eigenvecs;
	vector<double> vinitial(size, 0);
	for(int i = 0; i < size; i++) {
		eigenvecs.insert({i, vinitial});
	}

	// Create the identity matrix used as input for TQLI
	for(int i = 0; i < size; i++) {
		eigenvecs[i][i] = 1;
	}
	
	tqli(alpha, beta, size, eigenvecs);

	// result stores the eigenvalues of orginal Laplacian matrix, calculated by Matlab
	vector<double> result = {0, 1.381966, 2.381966, 3.61803, 4.61803};

	// Verify if the eigenvalues of the diagonalised the matrix are same as result
	for (int i = 0; i < size; i++) {
        if (abs(alpha[i] - result[i]) > 1e-5)
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

bool testPartition() {
	
	Graph g(100);
	Partition partition(g, 2, true);
	
	Analysis::cutEdgeVertexTable(g);	
	partition.printLapEigenvalues();
	
	return true;
}

bool testReothogonalisation() {

    Graph g(100);
    Partition partition1(g, 4, false);
	cout << "WITHOUT reorthogonalisation: " << endl;
	Analysis::cutEdgeVertexTable(g);
    partition1.printLapEigenvalues();
    partition1.printLapEigenMat();

    Partition partition2(g, 4, true);
	cout << "WITH reorthogonalisation: " << endl;
	Analysis::cutEdgeVertexTable(g);
    partition2.printLapEigenvalues();
    partition2.printLapEigenMat();

    return true;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testReadGraph
 *  Description:  Test the read function, can not verify the correctness in unit testing, but can test the performance
 * =====================================================================================
 */

bool testReadGraph() {
	Graph g;
	ifstream In("test_read_20.dot");
	if (!In.is_open()) {
		std::cerr << "ERROR: Can't open the file" << endl;
		exit(-1);
	}
	g.readDotFormat(In);
	
	g.printDotFormat();
	Partition partition(g, 4, true);
	g.printDotFormat();
	
	return true;
}

bool testReadGraphWithColour() {
	Graph g;
	ifstream In("test_read_20.dot");
	if (!In.is_open()) {
		std::cerr << "ERROR: Can't open the file" << endl;
		exit(-1);
	}
	g.readDotFormatWithColour(In);	

	g.printDotFormat();

	return true;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testCutEdgePercent
 *  Description:  Test the percentage of edges have been cut
 * =====================================================================================
 */

bool testCutEdgePercent() {
	Graph g;
	ifstream In("test_read_20.dot");
	if (!In.is_open()) {
		std::cerr << "ERROR: Can't open the file" << endl;
		exit(-1);
	}
	g.readDotFormatWithColour(In);	

	double cut_edge_percent = Analysis::cutEdgePercent(g);
	//cout << "cut_edge_percent = " << cut_edge_percent << endl;
	
	if (std::abs(cut_edge_percent - 0.142857) > 1e-5) {
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

bool testCutEdgeVertexTable() {
	Graph g;
	ifstream In("test_table_1000.dot");
	if (!In.is_open()) {
		std::cerr << "ERROR: Can't open the file" << endl;
		exit(-1);
	}
	g.readDotFormatWithColour(In);
	Analysis::cutEdgeVertexTable(g);
	
	/*-----------------------------------------------------------------------------
 	* Basic info of the graph
	*-----------------------------------------------------------------------------
	Vertices:   20
	Edges:      36
	Subgraphs:  4
	*-----------------------------------------------------------------------------
 	* Number of nodes in each subgraph
	*-----------------------------------------------------------------------------
	Colour:		0	1	2	3
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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testManuallyPartition
 *  Description:  Test the function manuallyPartition in analysis.cc
 * =====================================================================================
 */

bool testManuallyPartition() {
	Graph g;
	ifstream In("test_table_1000.dot");
	if (!In.is_open()) {
		std::cerr << "ERROR: Can't open the file" << endl;
		exit(-1);
	}
	g.readDotFormatWithColour(In);	
	//g.printDotFormat();
	//Partition partition(g, 4, true);
	Analysis::cutEdgeVertexTable(g);

	Analysis::manuallyPartition(g);
	g.outputDotFormat("Output_1000_manuallyPartition.dot");
	Analysis::cutEdgeVertexTable(g);

	return true;
}

}

int main() {
	//cout << Tests::testTqli() << endl;;
	//cout << Tests::testPartition() << endl;
	//cout << Tests::testLanczos() << endl;
	//cout << Tests::testCutEdgePercent() << endl;
	//Tests::testReadGraph();
	//Tests::testReadGraphWithColour();
	//Tests::testManuallyPartition();
	//Tests::testPartition();
	//Tests::testCutEdgeVertexTable();
	Tests::testReothogonalisation();

	return 0;
}
