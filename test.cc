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
	vec[0] = 1;
	//vec[0] = 0.7;
	//vec[1] = -2.3;
	//vec[2] = 15.6;
	//vec[3] = 21.6;
	//vec[4] = -7.2;

	//vec[0] = 0.5;
	//vec[1] = 0.5;
	//vec[2] = 0;
	//vec[3] = 0.5;
	//vec[4] = 0.5;

	// Normalise the initial vector
	cout << "Norm of vec = " << norm(vec) << endl;
	for (unsigned int i = 0; i < G.size(); i++)
		//vec[i] = vec[i]/norm(vec);
		vec[i] = vec[i]*(1/norm(vec));

	cout << "Norm of vec = " << norm(vec) << endl;
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
	int size = G.size();

	
	map<int, vector<double>> lanczos_vecs;

	map<pair<int,int>, double> trimat = constructTriMat(G, vec, alpha, beta, lanczos_vecs);
	beta.push_back(0);

	cout << "vector alpha: " << endl;
	for (const double& x:alpha)
		cout << x << " ";
	cout << endl;

	cout << "vector beta: " << endl;
	for (const double& x:beta)
		cout << x << " ";
	cout << endl;

	
	cout << "triangular matrix: " << endl;
	cout << "sizeoftrimat: " << trimat.size() << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << trimat[make_pair(i,j)] << "\t";
		cout << endl;
	}
    
	// TQLI test
	map<pair<int,int>, double> tri_eigen_vecs;
	for(int i = 0; i < size; i++) {
		tri_eigen_vecs[make_pair(i,i)] = 1;
	}

	//cout << "eigenvector matrix: " << endl;
	//cout << "sizeofeigenvec: " << tri_eigen_vecs.size() << endl;
	//for (int i = 0; i < size; i++) {
	//	for (int j = 0; j < size; j++)
	//		cout << tri_eigen_vecs[make_pair(i,j)] << "\t";
	//	cout << endl;
	//}
	//cout << endl;

	tqli(alpha, beta, size, tri_eigen_vecs);

	cout << "eigenvector matrix: " << endl;
	cout << "sizeofeigenvec: " << tri_eigen_vecs.size() << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << tri_eigen_vecs[make_pair(i,j)] << "\t";
		cout << endl;
	}

	cout << "vector alpha(eigenvalues): " << endl;
	for (const double& x:alpha)
		cout << x << " ";
	cout << endl;


  	//Calculate the eigenvectors of Laplacian matrix

	// The Lanczos vector

	vector<double> v0(G.size(), 0);
	v0[0] = 1;
	cout << "lanczos vecs: " << endl;
	for (int i = 0; i< size; i++) {
		auto it = lanczos_vecs.find(i);
		for (const double& x:it->second)
			cout << x << " ";
		cout << endl;
	}


	/*-----------------------------------------------------------------------------
	 *  Calculate the second eigenvector of the original matrix
	 *-----------------------------------------------------------------------------*/

/*
	vector<double> tri_second_vec(size, 0);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (j == 1)
				tri_second_vec[i] = tri_eigen_vecs[make_pair(i,j)];
		}
	}
	cout << "first eigenvec: " << endl;
	for (const double& x:tri_second_vec) {
		cout << x <<  " ";	
	}
	cout << endl;

	vector<double> original_second_vec(size, 0);
	for	(int i = 0; i < size; i++) {
		for	(int j = 0; j < size; j++) {
			original_second_vec[i] += lanczos_vecs[j][i] * tri_second_vec[j];
		}
	}

	cout << "Second original eigenvec: " << endl;
	for (const double& x:original_second_vec) {
		cout << x/original_second_vec[size-1] <<  " ";	
	}
	cout << endl;
*/

	/*-----------------------------------------------------------------------------
	 *  Calculate all the eigenvectors of the original matrix
	 *-----------------------------------------------------------------------------*/

	// Lanczos vector * Trivectors
	map<pair<int,int>, double> lapacian_vecs;
	for (int k = 0; k < size; k++) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
			lapacian_vecs[make_pair(k,i)] += lanczos_vecs[j][i] * tri_eigen_vecs[make_pair(j,k)];
			}
		}	
	}

	cout << "lapacian_vecs: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << lapacian_vecs[make_pair(i,j)]/lapacian_vecs[make_pair(i,size-1)] << " ";
		}
		cout << endl;
	}
	
	return 0;
}
