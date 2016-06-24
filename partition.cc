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

#include "partition.h"
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"

#define Sign(a) (a >= 0.0 ? 1:0)

using namespace std;

typedef std::vector<double> Vector;
typedef std::map<std::pair<int,int>, double> SparseMatrix;
typedef std::unordered_map<int, Vector> DenseMatrix;

void printSparseMatrix(SparseMatrix& sparse_matrix, int size);
void printDenseMatrix(DenseMatrix& dense_matrix);
void printVector(const Vector& vec);
//Vector getEigenVec(const Graph& g);
std::unordered_map<int, Vector> getEigenMatrix(const Graph& g);


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  utils
 *  Description:  Print vector, matrix for debug
 * =====================================================================================
 */

void printSparseMatrix(SparseMatrix& sparse_matrix, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << sparse_matrix[make_pair(i,j)] << "\t";
		}
		cout << endl;
	}
}

void printDenseMatrix(DenseMatrix& dense_matrix) {
	int row_size = dense_matrix.size();
	for (int row = 0; row < row_size; row++) {
		int col_size = dense_matrix[row].size();
		for (int col = 0; col < col_size; col++) {
			cout << dense_matrix[row][col] << " ";
		}
		cout << endl;
	}
}

void printVector(const Vector& vec) {
	for (const double& x:vec) {
		cout << x <<  " ";
	}
	cout << endl;
}

void printEigenvalues(const Graph& g) {
	int size = g.size();

	// Define the input arguments for Lanczos to construct tridiagonal matrix
	// Construct tridiagonal matrix using Lanczos algorithm

	//Lanczos<Vector> lanczos(g);
	Lanczos<std::vector<double>> lanczos(g);
	Vector alpha = lanczos.alpha;
	Vector beta = lanczos.beta;

	beta.push_back(0);

	// Define an identity matrix as the input for TQLI algorithm
	DenseMatrix tri_eigen_vecs;
	Vector vinitial(size,0);
	for(int i = 0; i < size; i++)	tri_eigen_vecs[i] = vinitial;
	for(int i = 0; i < size; i++)	tri_eigen_vecs[i][i] = 1;

	// Calculate the eigenvalues and eigenvectors of the tridiagonal matrix
	//tqli(alpha, beta, size, tri_eigen_vecs);
	tqli(alpha, beta, size, tri_eigen_vecs);

	ofstream Output;
	Output.open("eigenvalues.dat");

	for (int i = 0; i < size; i++)	{
		Output << i << " " << alpha[i] << endl;
	}

	Output.close();
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getEigenVec
 *  Description:  Return the fiedler vector using lanczos and tqli functions
 * =====================================================================================
 */

Vector getEigenVec(DenseMatrix& lanczos_vecs, DenseMatrix& tri_eigen_vecs, const int& vector_index, const int& size) {
	// Find the eigenvector from the vector matrix of the Tridiagonal matrix
	Vector tri_eigen_vec(size, 0);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (j == vector_index) {
				tri_eigen_vec[i] = tri_eigen_vecs[i][j];
			}
		}
	}

	// Calculate the corresponding Laplacian vector
	Vector laplacian_vector(size, 0);
	for	(int i = 0; i < size; i++) {
		for	(int j = 0; j < size; j++) {
			laplacian_vector[i] += lanczos_vecs[j][i] * tri_eigen_vec[j];
		}
	}

	return laplacian_vector;
} 


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getEigenMatrix
 *  Description:  Get all the eigenvectors of the original matrix.
 * =====================================================================================
 */

DenseMatrix getEigenMatrix(const Graph& g, Vector& eigenvalues) {

	int size = g.size();

	// Define the input arguments for Lanczos to construct tridiagonal matrix

	// Construct tridiagonal matrix using Lanczos algorithm
	Lanczos<Vector> lanczos(g);
	Vector alpha = lanczos.alpha;
	Vector beta = lanczos.beta;

	beta.push_back(0);

	// Define an identity matrix as the input for TQLI algorithm
	DenseMatrix tri_eigen_vecs;
	DenseMatrix laplacian_eigen_vecs;

	Vector vinitial(size,0);
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
			laplacian_eigen_vecs[k][i] += lanczos.lanczos_vecs[j][i] * tri_eigen_vecs[j][k];
			}
		}	
	}

	return laplacian_eigen_vecs;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  partition
 *  Description:  Partition the graph into multiple subgraphs by only calculating the corresponding laplacian eigenvectors
 * =====================================================================================
 */

void partition(const Graph& g, const int subgraphs) {

	int size = g.size();
	int num_of_eigenvec = log2(subgraphs);

	// Construct tridiagonal matrix using Lanczos algorithm
	Lanczos<Vector> lanczos(g);
	Vector alpha = lanczos.alpha;
	Vector beta = lanczos.beta;

	beta.push_back(0);

#ifdef Debug
	cout << endl;
	cout << "triangular matrix: " << endl;
	printSparseMatrix(trimat, size);
#endif

	// Define an identity matrix as the input for TQLI algorithm
	DenseMatrix tri_eigen_vecs;
	Vector vinitial(size,0);
	for(int i = 0; i < size; i++)	tri_eigen_vecs[i] = vinitial;
	for(int i = 0; i < size; i++)	tri_eigen_vecs[i][i] = 1;

	// Calculate the eigenvalues and eigenvectors of the tridiagonal matrix
	tqli(alpha, beta, size, tri_eigen_vecs);

	// Find the index of the nth smallest eigenvalue (fiedler vector) of the eigenvalues vector "alpha" 
	int vector_index = 0;
	DenseMatrix partial_laplacian_eigen_vecs;

	unordered_multimap<double, int> hashmap;
	for (unsigned int i = 0; i < alpha.size(); i++)
		hashmap.insert({alpha[i], i});

	Vector auxiliary_vec = alpha;
	sort(auxiliary_vec.begin(), auxiliary_vec.end());
	for (int i = 0; i < num_of_eigenvec; i++) {
		auto it = hashmap.find(auxiliary_vec[i+1]);
		vector_index = it->second;
		hashmap.erase(it);
		partial_laplacian_eigen_vecs[i] = getEigenVec(lanczos.lanczos_vecs, tri_eigen_vecs, vector_index, size);
	}

#ifdef Debug
	printDenseMatrix(partial_laplacian_eigen_vecs);
#endif

	for (int vertex = 0; vertex < size; vertex++) {
		int colour = 0;
		for (int row = 0; row < num_of_eigenvec; row++) {
			colour += pow(2, row) * Sign(partial_laplacian_eigen_vecs[row][vertex]);
		}
		g.setColour(vertex, colour);
	}
}

/*-----------------------------------------------------------------------------
 *  Modified partition function to for multiple partitioning by calculating all laplacian eigenvectors
 *-----------------------------------------------------------------------------*/

void multiPartition(const Graph& g) {
	int size = g.size();

	Vector eigenvalues;
	DenseMatrix laplacian_eigen_vecs = getEigenMatrix(g, eigenvalues);

#ifdef Debug
	// Print all the eigenvalues of the tridiagonal/laplacian matrix
	cout << "laplacian eigenvalues: " << endl;
	printVector(eigenvalues);

	// Print all the eigenvectors in row
	cout << "laplacian_eigen_vecs: " << endl;
	printDenseMatrix(laplacian_eigen_vecs);
#endif

	// eigenvalues_index_sort stores the original index of the sorted eigenvalues 
	vector<int> eigenvalues_index_sort;
	unordered_multimap<double, int> hashmap;
	for (int i = 0; i < size; i++)
		hashmap.insert({eigenvalues[i], i});

	Vector auxiliary_vec = eigenvalues;
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
			cout << Sign(laplacian_eigen_vecs[row][col])<< " ";
		}
		cout << endl;
	}
#endif

	int num = 3; // number of eigenvectors to partition the graph, start from the second smallest.
	for (int vertex = 0; vertex < size; vertex++) {
		int colour = 0;
		for (int eigenvec_index = 1; eigenvec_index <= num; eigenvec_index++) {
			int row = eigenvalues_index_sort[eigenvec_index];
			colour += pow(2, eigenvec_index - 1) * Sign(laplacian_eigen_vecs[row][vertex]);
		}
		g.setColour(vertex, colour);
	}
}
