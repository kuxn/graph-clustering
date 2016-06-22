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
#include <map>
#include <unordered_map>
#include "graph.h"

#define Sign(a) (a >= 0.0 ? 1:0)

void printSparseMatrix(map<pair<int,int>, double>& sparse_matrix, int size);
void printDenseMatrix(unordered_map<int, vector<double>>& dense_matrix);
void printVector(const vector<double>& vec);
vector<double> getEigenVec(const Graph& g);
unordered_map<int, vector<double>> getEigenMatrix(const Graph& g);
void partition(const Graph& g);
void multiPartition(const Graph& g);

#endif
