/*
 * =====================================================================================
 *
 *       Filename:  lanczos.h
 *
 *    Description:  The interface of lanczos algorithm
 *        Created:  06/29/2016 21:22:58
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef LANCZOS_H_
#define LANCZOS_H_

#include <vector>
#include <unordered_map>
#include <boost/mpi.hpp>
#include "graph.h"

template<typename Vector, typename T>
class Lanczos {
    private:
        boost::mpi::communicator world;
        Vector init(const Graph& g);
        Vector multGraphVec(const Graph& g, const Vector& vec);
        inline T dot(const Vector& v1, const Vector& v2);
        inline T norm(const Vector& vec);
        inline void gramSchmidt(const int& iter, Vector& v);
        const int getIteration(const int& num_of_eigenvec, const int & global_size);
    
        std::unordered_map<int, std::vector<int>> halo_recv; // <rank, halo_neighbours to receive>
        std::unordered_map<int, std::vector<int>> halo_send; // <rank, halo_neighbours to send>
        void haloInit(const Graph& g);
        void haloUpdate(const Graph& g, Vector& v_local, Vector& v_halo);
    
    public:
        Lanczos(const Graph& g, const int& num_of_eigenvec, bool GramSchmidt);
    
        Vector alpha;
        Vector beta;
        std::vector<Vector> lanczos_vecs;
        void print_tri_mat();
};

#include "../src/lanczos.cc"
#endif
