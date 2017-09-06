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
private:
    static double cutEdgePercent(const Graph& g);

public:
    static void cutEdgeVertexTable(const Graph& g,
                                   const std::vector<double>& ritz_values);
    static void manuallyPartition(const Graph& g);
    static void outputTimes(const int& procs, const int& size,
                            const std::vector<double>& vec);
};

#endif
