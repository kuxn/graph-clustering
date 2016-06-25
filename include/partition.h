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

class Partition {

    private:
        typedef std::vector<double> Vector;
        typedef std::unordered_map<int, Vector> DenseMatrix;

        Vector laplacian_eigen_vec_;
        DenseMatrix laplacian_eigen_mat_;
        
        Vector getOneLapEigenVec(DenseMatrix& lanczos_vecs, DenseMatrix& tri_eigen_vecs, const int& vector_index);
        void getLapEigenMat(const Graph& g);

    public:
        Partition(const Graph& g, const int& subgraphs);
        void usingFullMat(const Graph& g, const int& subgraphs);

        void printLapEigenMat();
        void printLapEigenvalues();
        void outputLapEigenvalues();
};

#endif
