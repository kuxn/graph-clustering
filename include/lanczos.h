/*
 * =====================================================================================
 *
 *       Filename:  lanczos.h
 *
 *    Description:  Header file for lanczos algorithm
 *        Created:  06/11/2016 21:22:58
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

std::map<std::pair<int,int>, double> constructTriMat(const Graph& g, std::vector<double>& alpha, std::vector<double>& beta, std::unordered_map<int, std::vector<double>>& lanczos_vecs);

#endif
