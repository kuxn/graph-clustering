/*
 * =====================================================================================
 *
 *       Filename:  tqli.h
 *
 *    Description:  Header file for tqli.cc
 *		  Created:  06/14/2016 22:25:00
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef TQLI_H_
#define TQLI_H_

#include <cmath>
#include <vector>
#include <unordered_map>
#include <utility>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) // Magnitude of a times sign of b

using namespace std;

double SQR(double a); 
double pythag(double a, double b);
void tqli (vector<double>& d, vector<double>& e, int n, unordered_map<int, vector<double>>& z);

#endif



