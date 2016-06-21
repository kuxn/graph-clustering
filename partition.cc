/*
 * =====================================================================================
 *
 *       Filename:  partition.cc
 *
 *    Description:  Partition the graph according to the eigenvectors
 *		  Created:  06/17/2016 14:15:15
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <algorithm>
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"

using namespace std;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getEigenVec
 *  Description:  Return the second vector using lanczos and tqli functions
 * =====================================================================================
 */

vector<double> getEigenVec(const Graph& g) {

  int size = g.size();

	// Define the input arguments for Lanczos to construct tridiagonal matrix
	vector<double> alpha, beta;
	
	vector<double> v_initial(size, 0);
	for (int i = 0; i < size; i++) {
		v_initial[i] = drand48();
	}
	double normalise = norm(v_initial);
	for (int i = 0; i < size; i++) {
		v_initial[i] /= normalise;
	}

	//cout << "norm(v_initial) = " << norm(v_initial) << endl;

	unordered_map<int, vector<double>> lanczos_vecs;

	// Construct tridiagonal matrix using Lanczos algorithm
	map<pair<int,int>, double> trimat = constructTriMat(g, v_initial, alpha, beta, lanczos_vecs);
	beta.push_back(0);

#ifdef Debug
	cout << endl;
	cout << "triangular matrix: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << trimat[make_pair(i,j)] << "\t";
		cout << endl;
	}
#endif

	// Define an identity matrix as the input for TQLI algorithm
	unordered_map<int, vector<double>> tri_eigen_vecs;
	vector<double> vinitial(size,0);
	for(int i = 0; i < size; i++) {
		tri_eigen_vecs[i] = vinitial;
	}

	for(int i = 0; i < size; i++) {
		tri_eigen_vecs[i][i] = 1;
	}

	// Calculate the eigenvalues and eigenvectors of the tridiagonal matrix
	tqli(alpha, beta, size, tri_eigen_vecs);

	// Find the index of the second smallest eigenvalue of the eigenvalues vector "alpha" 
	int index_of_second_vec = 0;
	unordered_multimap<double, int> hashmap;
	for (unsigned int i = 0; i < alpha.size(); i++)
		hashmap.insert({alpha[i], i});

	vector<double> auxiliary_vec = alpha;
	sort(auxiliary_vec.begin(), auxiliary_vec.end());
	auto it = hashmap.find(auxiliary_vec[1]);
	index_of_second_vec = it->second;
	
#ifdef Debug
	cout << "index_of_second_vec = " << index_of_second_vec << endl;
#endif

	// Calculate the second eigenvector of original Laplacian matrix using 
	// the second eigenvector of the tridiagonal matrix computed by TQLI 
	vector<double> tri_second_vec(size, 0);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (j == index_of_second_vec)
				tri_second_vec[i] = tri_eigen_vecs[i][j];
		}
	}
	//cout << "Second tri eigenvec: " << endl;
	//for (const double& x:tri_second_vec) {
	//	cout << x <<  " ";	
	//}
	//cout << endl;
	
	vector<double> second_eigen_vec(size, 0);
	for	(int i = 0; i < size; i++) {
		for	(int j = 0; j < size; j++) {
			second_eigen_vec[i] += lanczos_vecs[j][i] * tri_second_vec[j];
		}
	}

#ifdef Debug
	cout << "Second eigen vector: " << endl;
	for (const double& x:second_eigen_vec) {
		cout << x <<  " ";	
	}
	cout << endl;
#endif

	return second_eigen_vec;
} 


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getEigenMatrix
 *  Description:  Get all the eigenvectors of the original matrix.
 * =====================================================================================
 */

unordered_map<int, vector<double>> getEigenMatrix(const Graph& g, vector<double>& eigenvalues) {

	int size = g.size();

	// Define the input arguments for Lanczos to construct tridiagonal matrix
	vector<double> alpha, beta;
	vector<double> v_initial(size, 0);
	for (int i = 0; i < size; i++) {
		v_initial[i] = drand48();
	}
	double normalise = norm(v_initial);
	for (int i = 0; i < size; i++) {
		//cout << "v_initial["<<i<<"] = " << v_initial[i] << endl;	
		v_initial[i] /= normalise;
	}

	unordered_map<int, vector<double>> lanczos_vecs;

	// Construct tridiagonal matrix using Lanczos algorithm
	map<pair<int,int>, double> trimat = constructTriMat(g, v_initial, alpha, beta, lanczos_vecs);
	beta.push_back(0);

#ifdef Debug
	cout << endl;
	cout << "vector alpha: " << endl;
	for (const double& x:alpha)
		cout << x << " ";
	cout << endl;

	cout << "vector beta: " << endl;
	for (const double& x:beta)
		cout << x << " ";
	cout << endl;

	cout << "triangular matrix: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << trimat[make_pair(i,j)] << "\t";
		cout << endl;
	}
#endif
	cout << "triangular matrix: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << trimat[make_pair(i,j)] << "\t";
		cout << endl;
	}

	// Define an identity matrix as the input for TQLI algorithm
	unordered_map<int, vector<double>> tri_eigen_vecs;
	unordered_map<int, vector<double>> laplacian_eigen_vecs;
	vector<double> vinitial(size,0);
	for(int i = 0; i < size; i++) {
		tri_eigen_vecs[i] = vinitial;	
		laplacian_eigen_vecs[i] = vinitial;
	}
	for(int i = 0; i < size; i++) {
		tri_eigen_vecs[i][i] = 1;
	}

	// Calculate the eigenvalues and eigenvectors of the tridiagonal matrix
	tqli(alpha, beta, size, tri_eigen_vecs);
	eigenvalues = alpha;

	// Calculate all the eigenvectors of original Laplacian matrix using 
	// the eigenvectors of the tridiagonal matrix computed by TQLI 
	for (int k = 0; k < size; k++) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
			laplacian_eigen_vecs[k][i] += lanczos_vecs[j][i] * tri_eigen_vecs[j][k];
			}
		}	
	}

#ifdef Debug
	// Print all the eigenvalues of the tridiagonal/laplacian matrix
	cout << "laplacian eigenvalues: " << endl;
	for (const double & x:alpha) {
		cout << x << " ";
	}
	cout << endl;

	// Print all the eigenvectors in row
	cout << "laplacian_eigen_vecs: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << laplacian_eigen_vecs[i][j] << " ";
		}
		cout << endl;
	}
#endif
	return laplacian_eigen_vecs;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  partition
 *  Description:  Partition the graph into two subgraphs according to the second eigenvector
 * =====================================================================================
 */
void partition(const Graph& g) {

    int numofvertex = g.size();
    vector<double> second_eigen_vector = getEigenVec(g);
	cout << "Undirected Graph {" << endl;
	for (int vertex = 0; vertex < numofvertex; vertex++) {
		cout << vertex << "[Partition="<<Sign(second_eigen_vector[vertex])<<"];" << endl;
	}

	for (int vertex = 0; vertex < numofvertex; vertex++) {
		auto it = g.find(vertex);
		for (const int& neighbour:it->second)
			cout << vertex << "--" << neighbour << " ;" << endl;
	}
	cout << "}" << endl;
}

/*-----------------------------------------------------------------------------
 *  Modified partition function to for multiple partitioning
 *-----------------------------------------------------------------------------*/

void multiPartition(const Graph& g) {
	int size = g.size();

	vector<double> eigenvalues;
	unordered_map<int, vector<double>> laplacian_eigen_vecs = getEigenMatrix(g, eigenvalues);

	// Print all the eigenvalues of the tridiagonal/laplacian matrix
	cout << "laplacian eigenvalues: " << endl;
	for (const double & x:eigenvalues) {
		cout << x << " ";
	}
	cout << endl;

	// Print all the eigenvectors in row
	cout << "laplacian_eigen_vecs: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << laplacian_eigen_vecs[i][j] << " ";
		}
		cout << endl;
	}

	// eigenvalues_index_sort stores the original index of the sorted eigenvalues 
	vector<int> eigenvalues_index_sort;
	unordered_multimap<double, int> hashmap;
	for (int i = 0; i < size; i++)
		hashmap.insert({eigenvalues[i], i});

	vector<double> auxiliary_vec = eigenvalues;
	sort(auxiliary_vec.begin(), auxiliary_vec.end());
	for (int i = 0; i < size; i++) {
		auto it = hashmap.find(auxiliary_vec[i]);
		eigenvalues_index_sort.push_back(it->second);
		hashmap.erase(it);
	}
	cout << endl;

	cout << "Original index of eigenvalues in sorted order: "; 
	for (const int& x:eigenvalues_index_sort)
		cout <<	x << " "; 
	cout << endl;

	cout << "eigenvalues in sorted order: " << endl; 
	for (const int& x:eigenvalues_index_sort)
		cout <<	eigenvalues[x] << " "; 
	cout << endl;

	cout << "Laplacian eigenvectors in sorted order: " << endl;
	for (const int& row:eigenvalues_index_sort) {
		for (int col = 0; col < size; col++) {
			cout << laplacian_eigen_vecs[row][col] << " ";
		}
		cout << endl;
	}

	cout << "Laplacian eigenvectors in sorted order: " << endl;
	for (const int& row:eigenvalues_index_sort) {
		for (int col = 0; col < size; col++) {
			//int index = laplacian_eigen_vecs[row][col] >= 0.0 ? 1:0;
			cout << Sign(laplacian_eigen_vecs[row][col])<< " ";
		}
		cout << endl;
	}

	cout << "Undirected Graph {" << endl;
	int num = 2; // number of eigenvectors to partition the graph, start from the second smallest.
	for (int vertex = 0; vertex < size; vertex++) {
		int index = 0;
		for (int eigenvec_index = 1; eigenvec_index <= num; eigenvec_index++) {
			int row = eigenvalues_index_sort[eigenvec_index];
			index += pow(2, eigenvec_index - 1) * Sign(laplacian_eigen_vecs[row][vertex]);
		}
		cout << vertex << "[Partition="<<index<<"];" << endl;
	}

	for (int vertex = 0; vertex < size; vertex++) {
		auto it = g.find(vertex);
		for (const int& neighbour:it->second)
			cout << vertex << "--" << neighbour << " ;" << endl;
	}
	cout << "}" << endl;
}

