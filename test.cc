/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  Main function
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

using namespace std;

int main() {

	Graph G;

	// Addedge function test
	G.addEdge(0,1);
	G.addEdge(0,4);
	G.addEdge(1,2);
	G.addEdge(2,3);
	G.addEdge(2,4);
	G.addEdge(3,4);
	
	//G.printDotFormat();
	G.printLaplacianMat();

	// Lanczos test
	vector<double> vec(G.size(), 0);

	//vec[0] = 0.7;
	//vec[1] = -2.3;
	//vec[2] = 15.6;
	//vec[3] = 21.6;
	//vec[4] = -7.2;
/*
	vec[0] = 0.5;
	vec[1] = 0.5;
	vec[2] = 0;
	vec[3] = 0.5;
	vec[4] = 0.5;
*/

	vec[0] = 1;

	cout << "Input vector: " << endl;
	for (const double& x:vec)
		cout << x << " ";
	cout << endl;

	vector<double> vec2 = multGraphVec(G, vec);

	cout << "Grpah * Vec: " << endl;
	for (const double& x:vec2)
		cout << x << " ";
	cout << endl;

	double dotprod = dot(vec, vec2);
	cout << "dotprod = " << dotprod << endl;
	cout << "norm fo vec2 = " << norm(vec2) << endl;  

	vector<double> alpha, beta;
	map<pair<int,int>, double> trimat = constructTriMat(G, vec, alpha, beta);
	beta.push_back(0);

	cout << "vector alpha: " << endl;
	for (const double& x:alpha)
		cout << x << " ";
	cout << endl;

	cout << "vector beta: " << endl;
	for (const double& x:beta)
		cout << x << " ";
	cout << endl;

	int size = G.size();
	
	cout << "triangular matrix: " << endl;
	cout << "sizeoftrimat: " << trimat.size() << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << trimat[make_pair(i,j)] << "\t";
		cout << endl;
	}
  
	// TQLI test

	map<pair<int,int>, double> eigenvec;
	for(int i = 0; i < size; i++) {
		eigenvec[make_pair(i,i)] = 1;
	}

	cout << "eigenvector matrix: " << endl;
	cout << "sizeofeigenvec: " << eigenvec.size() << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << eigenvec[make_pair(i,j)] << "\t";
		cout << endl;
	}
	cout << endl;

	vector<double> diagonal, subdiagonal;
	diagonal = {6, 8, 11, 13};
	subdiagonal = {8, 2, 1, 0};

	tqli(diagonal, subdiagonal, 4, eigenvec);

	cout << "vector diagonal(eigenvalues): " << endl;
	for (const double& x:diagonal)
		cout << x << " ";
	cout << endl;

	cout << "vector subdiagonal(eigenvalues): " << endl;
	for (const double& x:subdiagonal)
		cout << x << " ";
	cout << endl;

	cout << "eigenvector matrix: " << endl;
	cout << "sizeofeigenvec: " << eigenvec.size() << endl;
	for (int i = 0; i < diagonal.size(); i++) {
		for (int j = 0; j < diagonal.size(); j++)
			cout << eigenvec[make_pair(i,j)] << "\t";
		cout << endl;
	}
	
	//tqli(alpha, beta, size, eigenvec);

	//cout << "eigenvector matrix: " << endl;
	//cout << "sizeofeigenvec: " << eigenvec.size() << endl;
	//for (int i = 0; i < size; i++) {
	//	for (int j = 0; j < size; j++)
	//		cout << eigenvec[make_pair(i,j)] << "\t";
	//	cout << endl;
	//}

	//cout << "vector alpha(eigenvalues): " << endl;
	//for (const double& x:alpha)
	//	cout << x << " ";
	//cout << endl;
	//cout << "vector beta: " << endl;
	//for (const double& x:beta)
	//	cout << x << " ";
	//cout << endl;

/*  Lanczos and TQLI test
	cout << "Matrix * Eigenvector test" << endl;
	vector<double> firstcol;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (j == 0)
			firstcol.push_back(eigenvec[make_pair(i,j)]);
		}
	}

	cout << "first col: " << endl;
	for (const double& x:firstcol)
		cout << x << " ";
	cout << endl;

	vector<double> result;
	result = multGraphVec(G, firstcol);

	cout << "orginal matrix * firstcol: " << endl;
	for (const double& x:result)
		cout << x << " ";
	cout << endl;
*/
	return 0;
}


