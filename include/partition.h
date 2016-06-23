/*
 * =====================================================================================
 *
 *       Filename:  partition.h
 *
 *    Description:  Header file for partition.cc
 *        Created:  06/17/2016 14:35:37
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include <vector>
#include <map>
#include "graph.h"

void printEigenvalues(const Graph& g);
void partition(const Graph& g, const int subgraphs);
void multiPartition(const Graph& g);

#endif
