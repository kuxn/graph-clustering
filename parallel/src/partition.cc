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

#include <boost/mpi/timer.hpp>

#include "partition.h"
#include "lanczos.h"
#include "tqli.h"
#include "vt_user.h"

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

Partition::Partition(const Graph& g, const int& subgraphs, bool GramSchmidt) {

    VT_TRACER("PARTITION");
    int size = g.globalSize();
    int num_of_eigenvec = log2(subgraphs);

    // Construct tridiagonal matrix using Lanczos algorithm
    boost::mpi::timer timer_lanczos;
    Lanczos<Vector, double> lanczos(g, GramSchmidt);

    if (g.rank() == 0) {
        cout << "In P0, Lanczos takes " << timer_lanczos.elapsed() << "s" << endl;
    }

    laplacian_eigenvalues_ = lanczos.alpha_global;
    Vector beta = lanczos.beta_global;
    beta.push_back(0.0);

#ifdef Debug
	cout << endl;
	if (g.rank() == 0) {
		cout << "lanczos matrix: " << endl;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				cout << lanczos.lanczos_vecs_global[i][j] << "\t";
			}
			cout << endl; 
		}
	}
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
    boost::mpi::timer timer_tqli;
    tqli(laplacian_eigenvalues_, beta, size, tri_eigen_vecs); 
    if (g.rank() == 0) {
        cout << "In P0, TQLI takes " << timer_tqli.elapsed() << "s" << endl;
    }

#ifdef Debug
	cout << "Eigenvalues after TQLI: ";
	for (auto x:laplacian_eigenvalues_) {
		cout << x << " ";
	}
	if (g.rank() == 0) {
		cout << "tridiagonal matrix: " << endl;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				cout << tri_eigen_vecs[i][j] << "\t";
			}
			cout << endl; 
		}
	}
#endif
	// Find the index of the nth smallest eigenvalue (fiedler vector) of the eigenvalues vector "alpha" 
    int vector_index = 0;

    unordered_multimap<double, int> hashmap;
    for (int i = 0; i < size; i++)
        hashmap.insert({laplacian_eigenvalues_[i], i});

    Vector auxiliary_vec = laplacian_eigenvalues_;
    sort(auxiliary_vec.begin(), auxiliary_vec.end());
    int fiedler_index = 1;
    for (int i = 0; i < num_of_eigenvec; i++) {
        auto it = hashmap.find(auxiliary_vec[fiedler_index]);
        while (abs(it->first) < 1e-2) {
            fiedler_index++;
            it = hashmap.find(auxiliary_vec[fiedler_index]);
        }
        fiedler_index++;
        vector_index = it->second;
        //if (g.rank() == 0) {
        //    cout << "eigenvalue used: " << it->first << ", Vector Index: " << vector_index <<endl;
        //}
        hashmap.erase(it);
        laplacian_eigen_mat_[i] = getOneLapEigenVec(lanczos.lanczos_vecs_global, tri_eigen_vecs, vector_index);
    }

    for (int vertex = 0; vertex < g.localSize(); vertex++) {
        int colour = 0;
        for (int row = 0; row < num_of_eigenvec; row++) {
            colour += pow(2, row) * Sign(laplacian_eigen_mat_[row][g.globalIndex(vertex)]);
        }
        g.setColour(g.globalIndex(vertex), colour);
    }
}

/*-----------------------------------------------------------------------------
 *  Modified partition function to for multiple partitioning by calculating all laplacian eigenvectors
 *-----------------------------------------------------------------------------*/

void Partition::usingFullMat(const Graph& g, const int& subgraphs, bool GramSchmidt) {

    int size = g.size();
    getLapEigenMat(g, GramSchmidt);

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
    cout << "Laplacian Eigenvalues: ";
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
    VT_TRACER("EigenVec");
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

void Partition::getLapEigenMat(const Graph& g, bool GramSchmidt) {

    int size = g.globalSize();

    // Construct tridiagonal matrix using Lanczos algorithm
    Lanczos<Vector, double> lanczos(g, GramSchmidt);
    laplacian_eigenvalues_ = lanczos.alpha_global;
    Vector beta = lanczos.beta_global;

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
                laplacian_eigen_mat_[k][i] += lanczos.lanczos_vecs_global[j][i] * tri_eigen_vecs[j][k];
            }
        }	
    }
}
