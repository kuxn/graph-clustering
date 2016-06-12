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

using namespace std;

int main() {

	Graph G;

	G.addEdge(0,1);
	G.addEdge(0,4);
	G.addEdge(1,2);
	G.addEdge(2,3);
	G.addEdge(2,4);
	G.addEdge(3,4);
	G.addEdge(3,5);
	
	//G.printDotFormat();

	G.printLaplacianMat();

	// Lanczos test

	vector<double> vec(G.size(), 0);
	vec[2] = 1;

	for (const int& x:vec)
		cout << x << " ";

	cout << endl;


	vector<double> vec2 = multGraphVec(G, vec);

	for (const int& x:vec2)
		cout << x << " ";

	double dotprod = dot(vec, vec2);
	cout << "dotprod = " << dotprod << endl;
	
	cout << "norm fo vec2 = " << norm(vec2) << endl;  

	map<pair<int,int>, double> trimat = constructTriMat(G, vec);

	int size = G.size();
	
	cout << "triangular matrix: " << endl;
	cout << "sizeoftrimat: " << trimat.size() << endl;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << trimat[make_pair(i,j)] << "\t";
		cout << endl;
	}
  
	return 0;
}


