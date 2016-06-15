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
	int normret = 0;
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
	vector<double> vprime, v1, v2;
	v1 = v0; v2 = v0;
	map<pair<int, int>, double> trimat;
	int size = v0.size();
	int iter = 0;
	double alphaval = 0, betaval = 0;
	alpha.push_back(0);
	beta.push_back(0);
	beta.push_back(0);

	vprime = multGraphVec(G, v0);
	alphaval = dot(v0, vprime)/dot(v0, v0);
	alpha.push_back(alphaval);
	trimat[make_pair(iter, iter)] = alphaval;
	for (int index = 0; index < size; index++)
		v1[index] = vprime[index] - alphaval * v0[index];

	betaval = norm(v1); 
	beta.push_back(betaval);	
	trimat[make_pair(iter, iter+1)] = betaval;
	trimat[make_pair(iter+1, iter)] = betaval;
	cout << "v"<< iter <<"*v" << iter+1 << " = " << dot(v0, v1) << endl;
	iter++;

	while(iter < size) {
		vprime = multGraphVec(G, v1);
		alphaval = dot(v1, vprime)/dot(v1, v1);
		alpha.push_back(alphaval);
		trimat[make_pair(iter, iter)] = alphaval;

		for (int index = 0; index < size; index++)
			v2[index] = vprime[index] - alphaval * v1[index] - betaval * v0[index];

		betaval = norm(v2); 
		beta.push_back(betaval);	
		trimat[make_pair(iter, iter+1)] = betaval;
		trimat[make_pair(iter+1, iter)] = betaval;
		v0 = v1; v1 = v2;
		cout << "v"<< iter <<"*v" << iter+1 << " = " << dot(v0, v1) << endl;
		iter++;
		cout << endl;
	}
	beta.erase(beta.end()-1); // Remove the last element
	trimat.erase(make_pair(iter-1, iter));
	trimat.erase(make_pair(iter, iter-1));
	return trimat;
}	
