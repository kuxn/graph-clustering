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

#include <iostream> 
#include <stdexcept>
#include <utility> 
#include <cmath> 
#include "lanczos.h"

using namespace std;

vector<double> multGraphVec(const Graph& g, const vector<double>& vec);
double dot(const vector<double>& v1, const vector<double>& v2);
double norm(const vector<double>& vec);



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  multGraphVec
 *  Description:  The first component of Lanczos iteration fomular, Laplacian matrix * vector
 * =====================================================================================
 */

vector<double> multGraphVec(const Graph& g, const vector<double>& vec) {
	vector<double> prod;
	if (g.size() != vec.size()) throw std::length_error("The sizes don't match.");
	int numofvertex = vec.size();
	for (int vertex = 0; vertex < numofvertex; vertex++) {
		auto it = g.find(vertex);
		double temp = 0;
		if (!it->second.empty())
			for (const int& neighbour:it->second)
				temp += vec[neighbour];

		prod.push_back(it->second.size() * vec[vertex] - temp);
	}
	return prod;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dot
 *  Description:  Dot product for the calcualtion of alpha for Lanczos
 * =====================================================================================
 */

double dot(const vector<double>& v1, const vector<double>& v2) {
	if (v1.size() != v2.size())	throw std::length_error("The vector sizes don't match.");
	int size = v1.size();
	double dotprod = 0;
	for (int index = 0; index < size; index++) {
		dotprod += v1[index] * v2[index];	
	}
	return dotprod;
}

double norm(const vector<double>& vec) {
	double normret = 0;
	for (const double& value:vec) {
		normret += value * value;
	}
	return sqrt(normret);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  constructTriMat
 *  Description:  The triangular matrix calculated by Lanczos
 * =====================================================================================
 */

/*-----------------------------------------------------------------------------
 *  Signiture of the funtion:
 * 	input:
 *		G contains all the elements of the original graph
 * 		v0[0..n-1] is the initial arbitrary vector whose norm is 1
 *		alpha[0..n-1] is an empty vector
 *		beta[0..n-1] is an empty vector
 *		lanczos_vecs[0..n-1][0..n-1] is an empty matrix
 *
 *	output:
 *		alpha[0..n-1] returns all the diagonal elements of the tridiagonal matrix (diagonalised)
 *		beta[0..n-1] returns all the subdiagonal elements of the tridiagonal matrix
 *		trimat[0..n-1][0..n-1] returns the tridiagonal matrix
 *		lanczos_vecs[0..n-1][0..n-1] the kth row returns the kth Lanczos vector
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 *  Modified Lanczos algorithm with reorthogonalisation by Gram–Schmidt
 *-----------------------------------------------------------------------------*/

map<pair<int,int>, double> constructTriMat(const Graph& g, vector<double>& alpha, vector<double>& beta, unordered_map<int, vector<double>>& lanczos_vecs) {


    int size = g.size();

	vector<double> v0(size, 0);
	for (int i = 0; i < size; i++)  v0[i] = drand48();

	double normalise = norm(v0);
	for (int i = 0; i < size; i++)  v0[i] /= normalise;

	vector<double> t = v0;
	vector<double> v1 = v0;
    vector<double> w;
   	map<pair<int, int>, double> trimat;

	double alpha_val = 0, beta_val = 0;
	lanczos_vecs[0] = v0;

	for (int iter = 1; iter < size; iter++) {
        w = multGraphVec(g, v1);
		alpha_val = dot(v1, w);
		//cout << "dot(v1, v1) = " << dot(v1, v1) << endl;
		alpha.push_back(alpha_val);
		trimat[make_pair(iter-1, iter-1)] = alpha_val;

		for (int index = 0; index < size; index++)
			t[index] = w[index] - alpha_val * v1[index] - beta_val * v0[index];

		beta_val = norm(t); 
		if (abs(beta_val) < 1e-5) 
		try { throw std::runtime_error("Value of beta is close to 0"); }
		catch (std::runtime_error& e) { std::cerr << "ERROR: " << e.what() << endl; }
		//cout << "beta_val: " << beta_val << endl;

		beta.push_back(beta_val);	
		trimat[make_pair(iter-1, iter)] = beta_val;
		trimat[make_pair(iter, iter-1)] = beta_val;

		//v0 = v1;
		for (int index = 0; index < size; index++)
			v1[index] = t[index]/beta_val;

		lanczos_vecs[iter] = v1;

		for (int k = 1; k <= iter; k++) {
			//cout << "i - norm of lanczos_vecs["<<k<<"] = " << norm(lanczos_vecs[k]) << endl;
			for (int i = 0; i < k; i++) {
				double reorthog_dot_product = dot(lanczos_vecs[i], lanczos_vecs[k]);
				for (int j = 0; j < size; j++) {
					lanczos_vecs[k][j] -= reorthog_dot_product * lanczos_vecs[i][j];
				}
			}
			double normalise = norm(lanczos_vecs[k]);
			for (int j = 0; j < size; j++) lanczos_vecs[k][j] /= normalise;
			//cout << "norm of lanczos_vecs["<<iter<<"] = " << norm(lanczos_vecs[iter]) << endl;
			for (int i = 0; i < k; i++) {
			//cout << "DOT: " << i << k << " "<< dot(lanczos_vecs[i], lanczos_vecs[k]) << endl;
			}
		}
		v0 = lanczos_vecs[iter-1];
		v1 = lanczos_vecs[iter];
		
		//Verify the dot product of v0 and v1 which is supposed to be 0
        double dot_product = dot(v0, v1);
#ifdef Debug		
		cout << "v"<< iter-1 <<"*v" << iter << " = " << dot_product << endl;
		cout << endl;
#endif
        if (abs(dot_product) > 1e-5) 
		try { throw std::runtime_error("Need reorthogonalise"); }
		catch (std::runtime_error& e) { std::cerr << "ERROR: " << e.what() << endl; }
	}
	w = multGraphVec(g, v1);
	alpha_val = dot(v1, w);
	alpha.push_back(alpha_val);
	trimat[make_pair(size-1, size-1)] = alpha_val;
	cout << "Lanczos algorithm is done." << endl;
	return trimat;
}
