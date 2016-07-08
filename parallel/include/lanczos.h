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
#include <map>
#include <boost/mpi.hpp>
#include "graph.h"

template<typename Vector, typename T>
class Lanczos {
    private:
        Vector& init(Vector& v_local, const Graph& g);
        Vector multGraphVec(const Graph& g, const Vector& vec);
        inline T dot(boost::mpi::communicator& world, const Vector& v1, const Vector& v2);
        inline T dot_local(const Vector& v1, const Vector& v2);
        inline T norm(const Vector& vec);
        inline Vector& normalise(Vector& vec);
        inline void gramSchmidt(int& iter, const int& size);
        
        std::unordered_map<int, std::unordered_set<int>> halo_recv; // <rank, halo_neighbours to receive>
        std::unordered_map<int, std::unordered_set<int>> halo_send; // <rank, halo_neighbours to send>
        void haloInit(boost::mpi::communicator& world, const Graph& g);
        void haloUpdate(boost::mpi::communicator& world, const Graph& g, Vector& v_local, Vector& v_halo);

    public:
        Lanczos(boost::mpi::communicator& world, const Graph& g, bool reorthogonalisation);

        Vector alpha_global;
        Vector beta_global;
        std::unordered_map<int, Vector> lanczos_vecs_global;
        std::unordered_map<int, Vector> lanczos_vecs_local;
        void print_tri_mat();

        void transform(boost::mpi::communicator& world, const Graph& g);
};

#include "../src/lanczos.cc"
#endif
