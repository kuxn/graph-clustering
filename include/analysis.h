/*
 * =====================================================================================
 *
 *       Filename:  analysis.h
 *
 *    Description:  Header file for analysis.cc
 *        Created:  06/22/2016 11:49:29
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "graph.h"

double cutEdgePercent(const Graph& g);
void cutEdgeNodeTable(const Graph& g, const int& subgraphs);
void manuallyPartition(const Graph& g, const int& subgraphs);

#endif   
