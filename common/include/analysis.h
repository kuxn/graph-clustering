/**
 * @file analysis.h
 * @brief Header file for analysis.cc
 * @author Ken Hu, xnchnhu@gmail.com
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <vector>
#include "graph.h"

class Analysis
{
public:
    static double cutEdgePercent(const Graph& g);
    static void cutEdgeVertexTable(const Graph& g,
                                   const std::vector<double>& ritzValues);
    static void manuallyPartition(const Graph& g);
    static void outputTimes(const int& procs, const int& size,
                            const std::vector<double>& vec);
    // For serial
    static void randomPartition(const Graph& g, const int& colours);
    static void evenPartition(const Graph& g, const int& colours);
    static void outputTimes(const int& numOfVertices, const std::vector<double>& vec);
    //static void benchmarks(bool enableGramSchmidt);
};

#endif
