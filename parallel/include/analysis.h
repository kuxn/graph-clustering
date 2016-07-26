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
#include <vector>

class Analysis {
    
    private:
        static double cutEdgePercent(const Graph& g);
    public:
        static void cutEdgeVertexTable(const Graph& g, const std::vector<double>& ritz_values);
        static void manuallyPartition(const Graph& g);
        static void outputTimes(const int& procs, const std::vector<double>& vec);
};

#endif   
