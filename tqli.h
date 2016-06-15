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
#include <map>
#include <utility>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) // Magnitude of a times sign of b

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) // Square a double value

using namespace std;

double pythag(double a, double b);
void tqli (vector<double>& d, vector<double>& e, int n, map<pair<int, int>, double>& z);

#endif



