/*
 * =====================================================================================
 *
 *       Filename:  partition.cc
 *
 *    Description:  Partition the graph according to the eigenvectors
 *        Created:  06/17/2016 14:15:15
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_map>

#include "partition.h"
#include "lanczos.h"
#include "tqli.h"

#define Sign(a) (a >= 0.0 ? 1:0)

using namespace std;

typedef std::vector<double> Vector;
typedef std::unordered_map<int, Vector> DenseMatrix;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Constructor
 *  Description:  Partition the graph into multiple subgraphs by only calculating the corresponding laplacian eigenvectors
 * =====================================================================================
 */

Partition::Partition(const Graph& g, const int& subgraphs) {

	int size = g.size();
	int num_of_eigenvec = log2(subgraphs);

	// Construct tridiagonal matrix using Lanczos algorithm
	Lanczos<Vector, double> lanczos(g);
	laplacian_eigenvalues_ = lanczos.alpha;
	Vector beta = lanczos.beta;

	beta.push_back(0);

#ifdef Debug
	cout << endl;
	cout << "triangular matrix: " << endl;
	lanczos.print_tri_mat();
#endif

	// Define an identity matrix as the input for TQLI algorithm
	DenseMatrix tri_eigen_vecs;
	Vector vinitial(size,0);
	for(int i = 0; i < size; i++)	tri_eigen_vecs[i] = vinitial;
	for(int i = 0; i < size; i++)	tri_eigen_vecs[i][i] = 1;

	// Calculate the eigenvalues and eigenvectors of the tridiagonal matrix
	tqli(laplacian_eigenvalues_, beta, size, tri_eigen_vecs);

	// Find the index of the nth smallest eigenvalue (fiedler vector) of the eigenvalues vector "alpha" 
	int vector_index = 0;

	unordered_multimap<double, int> hashmap;
	for (int i = 0; i < size; i++)
		hashmap.insert({laplacian_eigenvalues_[i], i});

	Vector auxiliary_vec = laplacian_eigenvalues_;
	sort(auxiliary_vec.begin(), auxiliary_vec.end());
	for (int i = 0; i < num_of_eigenvec; i++) {
		auto it = hashmap.find(auxiliary_vec[i+1]);
		vector_index = it->second;
		hashmap.erase(it);
		laplacian_eigen_mat_[i] = getOneLapEigenVec(lanczos.lanczos_vecs, tri_eigen_vecs, vector_index);
	}

#ifdef Debug
	printLapEigenMat();
#endif

	for (int vertex = 0; vertex < size; vertex++) {
		int colour = 0;
		for (int row = 0; row < num_of_eigenvec; row++) {
			colour += pow(2, row) * Sign(laplacian_eigen_mat_[row][vertex]);
		}
		g.setColour(vertex, colour);
	}
}

/*-----------------------------------------------------------------------------
 *  Modified partition function to for multiple partitioning by calculating all laplacian eigenvectors
 *-----------------------------------------------------------------------------*/

void Partition::usingFullMat(const Graph& g, const int& subgraphs) {

	int size = g.size();
	getLapEigenMat(g);

#ifdef Debug
	// Print all the eigenvalues of the tridiagonal/laplacian matrix
	cout << "laplacian eigenvalues: " << endl;
	printLapEigenvalues();

	// Print all the eigenvectors in row
	cout << "laplacian_eigen_mat_: " << endl;
	printLapEigenMat();
#endif

	// eigenvalues_index_sort stores the original index of the sorted eigenvalues 
	vector<int> eigenvalues_index_sort;
	unordered_multimap<double, int> hashmap;
	for (int i = 0; i < size; i++)
		hashmap.insert({laplacian_eigenvalues_[i], i});

	Vector auxiliary_vec = laplacian_eigenvalues_;
	sort(auxiliary_vec.begin(), auxiliary_vec.end());
	for (int i = 0; i < size; i++) {
		auto it = hashmap.find(auxiliary_vec[i]);
		eigenvalues_index_sort.push_back(it->second);
		hashmap.erase(it);
	}

#ifdef Debug
	cout << endl;
	cout << "Original index of eigenvalues in sorted order: "; 
	for (const int& x:eigenvalues_index_sort)
		cout <<	x << " "; 
	cout << endl;

	cout << "eigenvalues in sorted order: " << endl; 
	for (const int& x:eigenvalues_index_sort)
		cout <<	laplacian_eigenvalues_[x] << " "; 
	cout << endl;

	cout << "Laplacian eigenvectors in sorted order: " << endl;
	for (const int& row:eigenvalues_index_sort) {
		for (int col = 0; col < size; col++) {
			cout << laplacian_eigen_mat_[row][col] << " ";
		}
		cout << endl;
	}

	cout << "Laplacian eigenvectors in sorted order: " << endl;
	for (const int& row:eigenvalues_index_sort) {
		for (int col = 0; col < size; col++) {
			cout << Sign(laplacian_eigen_mat_[row][col])<< " ";
		}
		cout << endl;
	}
#endif

	int num = log2(subgraphs); // number of eigenvectors to partition the graph, start from the second smallest.
	for (int vertex = 0; vertex < size; vertex++) {
		int colour = 0;
		for (int eigenvec_index = 1; eigenvec_index <= num; eigenvec_index++) {
			int row = eigenvalues_index_sort[eigenvec_index];
			colour += pow(2, eigenvec_index - 1) * Sign(laplacian_eigen_mat_[row][vertex]);
		}
		g.setColour(vertex, colour);
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  utilities
 *  Description:  Print vector, matrix for debug
 * =====================================================================================
 */

void Partition::printLapEigenMat() {
	int row_size = laplacian_eigen_mat_.size();
	for (int row = 0; row < row_size; row++) {
		int col_size = laplacian_eigen_mat_[row].size();
		for (int col = 0; col < col_size; col++) {
			cout << laplacian_eigen_mat_[row][col] << " ";
		}
		cout << endl;
	}
}

void Partition::printLapEigenvalues() {
	for (const double& x:laplacian_eigenvalues_) {
		cout << x <<  " ";
	}
	cout << endl;
}

void Partition::outputLapEigenvalues() {

	string filename("eigenvalues_");
	filename += to_string(laplacian_eigenvalues_.size());
	filename += ".dat";

	ofstream Output(filename);

	for (unsigned int i = 0; i < laplacian_eigenvalues_.size(); i++)	{
		Output << i << " " << laplacian_eigenvalues_[i] << endl;
	}

	Output.close();
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getOneLapEigenVec
 *  Description:  Return the fiedler vector using lanczos and tqli functions
 * =====================================================================================
 */

Vector Partition::getOneLapEigenVec(DenseMatrix& lanczos_vecs, DenseMatrix& tri_eigen_vecs, const int& vector_index) {
	// Find the eigenvector from the vector matrix of the Tridiagonal eigenvector matrix (each column represents a vector)
	int size = tri_eigen_vecs.size(); // Size of vertices/vectors
	Vector tri_eigen_vec(size, 0);
	for (int row = 0; row < size; row++) {
		for (int col = 0; col < size; col++) {
			if (col == vector_index) {
				tri_eigen_vec[row] = tri_eigen_vecs[row][col];
			}
		}
	}

	// Calculate the corresponding Laplacian vector by Lanczos vectors(each row represents a vector, need transposing)
	Vector laplacian_vector(size, 0);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			laplacian_vector[i] += lanczos_vecs[j][i] * tri_eigen_vec[j];
		}
	}

	return laplacian_vector;
} 

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getLapEigenMat
 *  Description:  Get all the eigenvectors of the original matrix.
 * =====================================================================================
 */

void Partition::getLapEigenMat(const Graph& g) {

	int size = g.size();

	// Construct tridiagonal matrix using Lanczos algorithm
	Lanczos<Vector, double> lanczos(g);
	laplacian_eigenvalues_ = lanczos.alpha;
	Vector beta = lanczos.beta;

	beta.push_back(0);

	// Define an identity matrix as the input for TQLI algorithm
	DenseMatrix tri_eigen_vecs;

	Vector vinitial(size,0);
	for(int i = 0; i < size; i++) {
		tri_eigen_vecs[i] = vinitial;	
		laplacian_eigen_mat_[i] = vinitial;
	}
	for(int i = 0; i < size; i++) {
		tri_eigen_vecs[i][i] = 1;
	}

	// Calculate the eigenvalues and eigenvectors of the tridiagonal matrix
	tqli(laplacian_eigenvalues_, beta, size, tri_eigen_vecs);

	// Calculate all the eigenvectors of original Laplacian matrix using 
	// the eigenvectors of the tridiagonal matrix computed by TQLI 
	for (int k = 0; k < size; k++) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
			laplacian_eigen_mat_[k][i] += lanczos.lanczos_vecs[j][i] * tri_eigen_vecs[j][k];
			}
		}	
	}
}

