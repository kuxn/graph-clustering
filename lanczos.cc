/*
 * =====================================================================================
 *
 *       Filename:  lanczos.cpp
 *
 *    Description:  Lanczos algorithm
 *		  Created:  06/11/2016 15:33:49
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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  multGraphVec
 *  Description:  The first component of Lanczos iteration fomular
 * =====================================================================================
 */

vector<double> multGraphVec(const Graph& G, const vector<double>& vec) {
	vector<double> prod;
	if (G.size() != vec.size()) throw std::length_error("The sizes don't match.");
	int numofvertex = vec.size();
	for (int vertex = 0; vertex < numofvertex; vertex++) {
		auto it = G.find(vertex);
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

map<pair<int,int>, double> constructTriMat(const Graph& G, vector<double>& v0, vector<double>& alpha, vector<double>& beta) {
	vector<double> w, t, v1;
	t = v0; v1 = v0;
	map<pair<int, int>, double> trimat;
	int size = v0.size();
	int iter = 0;
	double alphaval = 0, betaval = 0;

	for (int iter = 1; iter < size; iter++) {
		w = multGraphVec(G, v1);
		alphaval = dot(v1, w);
		//cout << "dot(v1, v1) = " << dot(v1, v1) << endl;
		alpha.push_back(alphaval);
		trimat[make_pair(iter-1, iter-1)] = alphaval;

		for (int index = 0; index < size; index++)
			t[index] = w[index] - alphaval * v1[index] - betaval * v0[index];

		betaval = norm(t); 
		beta.push_back(betaval);	
		trimat[make_pair(iter-1, iter)] = betaval;
		trimat[make_pair(iter, iter-1)] = betaval;

		v0 = v1;

		for (int index = 0; index < size; index++)
			v1[index] = t[index]/betaval;

		cout << "v"<< iter <<"*v" << iter+1 << " = " << dot(v0, v1) << endl;
		cout << endl;
	}
	w = multGraphVec(G, v1);
	alphaval = dot(v1, w);
	alpha.push_back(alphaval);
	trimat[make_pair(size-1, size-1)] = alphaval;
	return trimat;
}	
