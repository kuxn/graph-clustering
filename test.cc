/*
 * =====================================================================================
 *
 *       Filename:  test.cc
 *
 *    Description:  Functions for unit tests
 *		  Created:  06/09/2016 18:12:52
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "test.h"
#include "partition.h"

namespace Tests {

using std::vector;
using std::unordered_map;
using std::map;
using std::cout;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  testTqli
 *  Description:  Test tqli function for the correctness of calculating eigenvalues, 
 * 				  correct eigenvalues come from Matlab.
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
        if (diagonal[i] - result[i] > 1e-5)
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
	// Initialise the input vector for lanczos algorithm
	vector<double> v_initial(size, 0);
	for (int i = 0; i < size; i++) {
		v_initial[i] = drand48();
	}
	double normalise = norm(v_initial);
	for (int i = 0; i < size; i++) {
		v_initial[i] /= normalise;
	}

	// Calculate the diagonal and subdiagonal vectors
	vector<double> alpha, beta;
	unordered_map<int, vector<double>> lanczos_vecs;
	map<pair<int,int>, double> trimat = constructTriMat(g, v_initial, alpha, beta, lanczos_vecs);
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
        if (alpha[i] - result[i] > 1e-5)
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
	
	Graph g;
	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(1,2);
	g.addEdge(2,3);
	g.addEdge(3,4);
	g.addEdge(3,5);
	g.addEdge(4,5);
	
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
	ifstream In;
	In.open("read_test.dot");

	if (In) g.readDotFormat(In);
    
    g.printDotFormat();

    return true;
}

bool testReadGraphWithColour() {
	Graph g;
	ifstream In;
	In.open("read_test.dot");

	if (In) g.readDotFormatWithColour(In);	

	g.printDotFormat();
	partition(g);

	return true;
}
}

int main() {
	//cout << Tests::testTqli() << endl;;
	//cout << Tests::testPartition() << endl;
	//Tests::testReadGraph();
	Tests::testReadGraphWithColour();
	return 0;

}
