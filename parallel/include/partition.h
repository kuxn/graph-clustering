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
#include <boost/mpi.hpp>
#include "graph.h"

class Partition {

    private:
        typedef std::vector<double> Vector;
        typedef std::unordered_map<int, Vector> DenseMatrix;

        Vector laplacian_eigenvalues_;
        DenseMatrix laplacian_eigen_mat_;
        
        Vector getOneLapEigenVec(DenseMatrix& lanczos_vecs, DenseMatrix& tri_eigen_vecs, const int& vector_index);
        void getLapEigenMat(boost::mpi::communicator& world, const Graph& g, bool reorthogonalisation);

    public:
        Partition() {}
        Partition(boost::mpi::communicator& world, const Graph& g, const int& subgraphs, bool reorthogonalisation);
        void usingFullMat(boost::mpi::communicator& world, const Graph& g, const int& subgraphs, bool reorthogonalisation);

        void printLapEigenMat();
        void printLapEigenvalues();
        void outputLapEigenvalues();
};

#endif
