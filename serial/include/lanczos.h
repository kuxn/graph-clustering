/*
 * =====================================================================================
 *
 *       Filename:  lanczos.h
 *
 *    Description:  The interface of lanczos algorithm
 *        Created:  06/11/2016 21:22:58
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef LANCZOS_H_
#define LANCZOS_H_

#include <vector>
#include <map>
#include "graph.h"

template<typename Vector, typename T>
class Lanczos {
    private:
        Vector init(const int& size);
        Vector multGraphVec(const Graph& g, const Vector& vec);
        const int getIteration(const int& num_of_eigenvec, const int & size);
        inline T dot(const Vector& v1, const Vector& v2);
        inline T norm(const Vector& vec);
        inline Vector& normalise(Vector& vec);
        inline void gramSchmidt(const int& iter, Vector& w);
        inline T l2norm(const Vector& alpha, const Vector& beta);

    public:
        Lanczos(const Graph& g, const int& num_of_eigenvec, bool GramSchmidt);

        Vector alpha;
        Vector beta;
        std::vector<Vector> lanczos_vecs;
        void print_tri_mat();
};

#include "../src/lanczos.cc"
#endif
