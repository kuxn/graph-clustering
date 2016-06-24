/*
 * =====================================================================================
 *
 *       Filename:  lanczos.h
 *
 *    Description:  Header file for lanczos algorithm
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

template<typename Vector>
class Lanczos {
    private:
        Vector& initialise(Vector& vec);
        Vector multGraphVec(const Graph& g, const Vector& vec);
        inline double dot(const Vector& v1, const Vector& v2);
        inline double norm(const Vector& vec);
        inline Vector& normalise(Vector& vec);
        inline void gramSchmidt(int& iter, int& size);

    public:
        Lanczos(const Graph& g);

        Vector alpha;
        Vector beta;
        std::unordered_map<int, Vector> lanczos_vecs;
        void print_tri_mat();
};

#include "../lanczos.cc"

#endif
