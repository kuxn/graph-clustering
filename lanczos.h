/*
 * =====================================================================================
 *
 *       Filename:  lanczos.h
 *
 *    Description:  Header file for lanczos algorithm
 *		  Created:  06/11/2016 21:22:58
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef LANCZOS_H_
#define LANCZOS_H_

#include <vector>
#include <map>
#include "graph.h"

using namespace std;

vector<double> multGraphVec(const Graph& G, const vector<double>& vec);
double dot(const vector<double>& v1, const vector<double>& v2);
double norm(const vector<double>& vec);
map<pair<int,int>, double> constructTriMat(const Graph& G, vector<double>& v0, vector<double>& alpha, vector<double>& beta);	

#endif
