/*
 * =====================================================================================
 *
 *       Filename:  lanczos.cpp
 *
 *    Description:  Lanczos algorithm
 *        Created:  06/11/2016 15:33:49
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef LANCZOS_CPP_
#define LANCZOS_CPP_

#include <iostream> 
#include <stdexcept>
#include <utility> 
#include <cmath> 
#include "lanczos.h"

using std::cout;
using std::endl;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Constructor
 *  Description:  The triangular matrix calculated by Lanczos
 * =====================================================================================
 */

/*-----------------------------------------------------------------------------
 *  Signiture of the funtion:
 * 	input:
 *		G contains all the elements of the original graph
 *
 *	output:
 *		alpha[0..n-1] returns all the diagonal elements of the tridiagonal matrix (diagonalised)
 *		beta[0..n-1] returns all the subdiagonal elements of the tridiagonal matrix
 *		tri_mat[0..n-1][0..n-1] returns the tridiagonal matrix
 *		lanczos_vecs[0..n-1][0..n-1] the kth row returns the kth Lanczos vector
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 *  Modified Lanczos algorithm with reorthogonalisation by Gram–Schmidt
 *-----------------------------------------------------------------------------*/

template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g, bool reorthogonalisation) {

    int size = g.size();
	Vector v0(size);
	v0 = initialise(v0);

	Vector t = v0, v1 = v0, w;

	T alpha_val = 0.0, beta_val = 0.0;
	lanczos_vecs[0] = v0;

	for (int iter = 1; iter < size; iter++) {
		w = multGraphVec(g, v1);
		T alpha_val = dot(v1, w);
		//cout << "dot(v1, v1) = " << dot(v1, v1) << endl;
		alpha.push_back(alpha_val);

		for (int index = 0; index < size; index++) {
			t[index] = w[index] - alpha_val * v1[index] - beta_val * v0[index];
		}

		beta_val = norm(t); 
		beta.push_back(beta_val);	
		if (std::abs(beta_val) < 1e-5) 
		try { throw std::runtime_error("Value of beta is close to 0: "); }
		catch (std::runtime_error& e) { 
			std::cerr << "ERROR: " << e.what(); 
			cout << "beta[" << iter-1 << "]: " << beta_val << endl;
		}

		if (!reorthogonalisation) {
			v0 = v1;
		}
		for (int index = 0; index < size; index++) {
			v1[index] = t[index]/beta_val;
		}

		lanczos_vecs[iter] = v1;
		if (reorthogonalisation) {
			gramSchmidt(iter, size);
			v0 = lanczos_vecs[iter-1];
			v1 = lanczos_vecs[iter];
		}
		
		//Verify the dot product of v0 and v1 which is supposed to be 0
		T dot_product = dot(v0, v1);
#ifdef Debug
		cout << "v"<< iter-1 <<"*v" << iter << " = " << dot_product << endl;
		cout << endl;
#endif
		if (std::abs(dot_product) > 1e-5) 
		try { throw std::runtime_error("Need reorthogonalise: "); }
		catch (std::runtime_error& e) { 
			std::cerr << "ERROR: " << e.what(); 
			cout << "v"<< iter-1 <<"*v" << iter << " = " << dot_product << endl;
		}
	}
	w = multGraphVec(g, v1);
	alpha_val = dot(v1, w);
	alpha.push_back(alpha_val);
	if (reorthogonalisation) {
		cout << "Lanczos algorithm WITH reorthogonalisation is done." << endl;
	} else {
		cout << "Lanczos algorithm WITHOUT reorthogonalisation is done." << endl; 
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  multGraphVec
 *  Description:  The first component of Lanczos iteration fomular, Laplacian matrix * vector
 * =====================================================================================
 */

template<typename Vector, typename T>
Vector Lanczos<Vector, T>::multGraphVec(const Graph& g, const Vector& vec) {
	Vector prod;
	if (g.size() != (int)vec.size())
	throw std::length_error("The sizes don't match.");

	int numofvertex = vec.size();
	for (int vertex = 0; vertex < numofvertex; vertex++) {
		auto it = g.find(vertex);
		T temp = 0.0;
		if (!it->second.empty())
		for (const int& neighbour:it->second)
			temp += vec[neighbour];

		prod.push_back(it->second.size() * vec[vertex] - temp);
	}
	return prod;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  gramSchmidt
 *  Description:  Reorthogonalisation
 * =====================================================================================
 */

template<typename Vector, typename T>
inline void Lanczos<Vector, T>::gramSchmidt(int& iter, int& size) {
	for (int k = 1; k <= iter; k++) {
		//cout << "i - norm of lanczos_vecs["<<k<<"] = " << norm(lanczos_vecs[k]) << endl;
		for (int i = 0; i < k; i++) {
			T reorthog_dot_product = dot(lanczos_vecs[i], lanczos_vecs[k]);
			for (int j = 0; j < size; j++) {
				lanczos_vecs[k][j] -= reorthog_dot_product * lanczos_vecs[i][j];
			}
		}
		//T normalise = norm(lanczos_vecs[k]);
		//for (int j = 0; j < size; j++) lanczos_vecs[k][j] /= normalise;
		normalise(lanczos_vecs[k]);
		//cout << "norm of lanczos_vecs["<<iter<<"] = " << norm(lanczos_vecs[iter]) << endl;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  utilities
 *  Description:  Vector operations
 * =====================================================================================
 */

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::dot(const Vector& v1, const Vector& v2) {
	if (v1.size() != v2.size())	
	throw std::length_error("The vector sizes don't match.");
	int size = v1.size();
	T dotprod = 0.0;
	for (int index = 0; index < size; index++) {
		dotprod += v1[index] * v2[index];	
	}
	return dotprod;
}

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::norm(const Vector& vec) {
	T normret = 0.0;
	for (const T& value:vec) {
		normret += value * value;
	}
	return sqrt(normret);
}

template<typename Vector, typename T>
inline Vector& Lanczos<Vector, T>::normalise(Vector& vec) {
	int size = vec.size();
	T normal = norm(vec);
	for (int i = 0; i < size; i++) {
		vec[i] /= normal;
	}
	return vec;
}

template<typename Vector, typename T>
Vector& Lanczos<Vector, T>::initialise(Vector& vec) {
	int size = vec.size();
	for (int i = 0; i < size; i++) {
		vec[i] = drand48();
	}
	T normalise = norm(vec);
	for (int i = 0; i < size; i++) {
		vec[i] /= normalise;
	}
	return vec;
}

template<typename Vector, typename T>
void Lanczos<Vector, T>::print_tri_mat() {
	int size = alpha.size();
	for (int row = 0; row < size; row++) {
		for (int col = 0; col < size; col++) {
			if (row == col)
				cout << alpha[row] << "\t";
			else if (col - row == 1)
				cout << beta[row] << "\t";
			else if (row - col == 1)
				cout << beta[col] << "\t";
			else
				cout << "0" << "\t";
		}
		cout << endl;
	}
}
#endif
