/*
 * =====================================================================================
 *
 *       Filename:  partition.h
 *
 *    Description:  Header file for partition.cc
 *		  Created:  06/17/2016 14:35:37
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include <vector>
#include "graph.h"

#define Sign(a) (a >= 0.0 ? 1:0)

vector<double> getEigenVec(const Graph& g);
unordered_map<int, vector<double>> getEigenMatrix(const Graph& g);
void partition(const Graph& g);
void multiPartition(const Graph& g);

#endif
