/*
 * =====================================================================================
 *
 *       Filename:  partition.h
 *
 *    Description:  The interface of partitioning of a graph
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

class Partition {

    private:
        typedef std::vector<double> Vector;
        typedef std::vector<Vector> DenseMatrix;

        Vector laplacian_eigenvalues_;
        DenseMatrix laplacian_eigen_mat_;

        Vector getOneLapEigenVec(DenseMatrix& lanczos_vecs, DenseMatrix& tri_eigen_vecs, const int& vector_index);
        inline int signMedian(double entry, double median);

    public:
        Partition() {}
        Partition(const Graph& g, const int& subgraphs, bool GramSchmidt);

        void printLapEigenMat();
        void printLapEigenvalues();
        void outputLapEigenvalues();
        std::vector<double> ritz_values;
        std::vector<double> times;
};

#endif
