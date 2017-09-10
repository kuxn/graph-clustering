/**
 * @file partition.h
 * @brief The interface of partitioning of a graph
 * @author Ken Hu, xnchnhu@gmail.com
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include <map>
#include <vector>
#include "graph.h"

class Partition
{
private:
    using DenseMatrix = std::vector<std::vector<double>>;

    std::vector<double> laplacianEigenvalues_;
    DenseMatrix laplacianEigenMatrix_;

    std::vector<double> getOneLapEigenVec(DenseMatrix& lanczosVectors,
                                          DenseMatrix& tridiagonalEigenVectors,
                                          const int& vectorIndex);
    inline int signMedian(double entry, double median);

public:
    Partition() {}
    Partition(const Graph& g, const int& subgraphs, bool GramSchmidt);

    void printLapEigenMat();
    void printLapEigenvalues();
    void outputLapEigenvalues();
    std::vector<double> ritzValues;
    std::vector<double> times;
};

#endif
