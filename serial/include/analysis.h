/**
 * @file analysis.h
 * @brief Header file for analysis.cc
 * @author Ken Hu, xnchnhu@gmail.com
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "graph.h"
#include <vector>

class Analysis {

    public:
        static double cutEdgePercent(const Graph& g);
        static void cutEdgeVertexTable(const Graph& g, const std::vector<double>& ritz_values);
        static void randomPartition(const Graph& g, const int& colours);
        static void evenPartition(const Graph& g, const int& colours);
        static void benchmarks(bool GramSchmidt);
        static void outputTimes(const int& size, const std::vector<double>& vec);
};

#endif   
