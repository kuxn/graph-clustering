/*
 * =====================================================================================
 *
 *       Filename:  lanczos.cc
 *
 *    Description:  Lanczos algorithm
 *        Created:  06/29/2016 21:33:49
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef LANCZOS_CPP_
#define LANCZOS_CPP_

#include <iostream>
#include <random>
#include <chrono>
#include <exception>
#include <utility>
#include <cmath>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/mpi/timer.hpp>

#include "lanczos.h"
//#include "vt_user.h"

namespace mpi = boost::mpi;
using std::cout;
using std::endl;

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Constructor
 *  Description:  The triangular matrix calculated by Lanczos
 * =====================================================================================
 */

/*-----------------------------------------------------------------------------
 *  Signiture of the funtion:
 * 	input:
 *		G contains all the elements of the original graph
 *
 *	output:
 *		alpha[0..n-1] returns all the diagonal elements of the tridiagonal matrix (diagonalised)
 *		beta[0..n-1] returns all the subdiagonal elements of the tridiagonal matrix
 *		tri_mat[0..n-1][0..n-1] returns the tridiagonal matrix
 *		lanczos_vecs_global[0..n-1][0..n-1] the kth row returns the kth Lanczos vector
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 *  Modified Lanczos algorithm with GramSchmidt by Gram–Schmidt
 *-----------------------------------------------------------------------------*/

#ifdef GS_
template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g_local, bool GramSchmidt) {

    //VT_TRACER("LANCZOS");
    int local_size = g_local.localSize();
    Vector v1_halo(g_local.globalSize());

    Vector v0_local = init(g_local);

    Vector v1_local = v0_local, w_local;
    T alpha_val_global = 0.0, beta_val_global = 0.0;

    lanczos_vecs_local[0] = v1_local;
    haloInit(g_local);

    int m = g_local.globalSize();
    //int m = 2 * std::sqrt(g_local.globalSize());

    for (int iter = 1; iter < m; iter++) {
        haloUpdate(g_local, v1_local, v1_halo);

        w_local = multGraphVec(g_local, v1_halo);

        alpha_val_global = dot(v1_local, w_local);
        alpha_global.push_back(alpha_val_global);

        for (int i = 0; i < local_size; i++) {
            w_local[i] = w_local[i] - alpha_val_global * v1_local[i] - beta_val_global * v0_local[i];
        }

        beta_val_global = sqrt(dot(w_local, w_local));
        beta_global.push_back(beta_val_global);

        if (std::abs(beta_val_global) < 1e-5)
            try {
                throw std::runtime_error("Value of beta is close to 0: ");
            }
            catch (std::runtime_error& e) {
                std::cerr << "ERROR: " << e.what();
                cout << "beta[" << iter-1 << "]: " << beta_val_global << endl;
            }

        for (int i = 0; i < local_size; i++) {
            v1_local[i] = w_local[i]/beta_val_global;
        }

        if (GramSchmidt) {
            gramSchmidt(iter, v1_local);
        }

        lanczos_vecs_local[iter] = v1_local;
        v0_local = lanczos_vecs_local[iter-1];

        //Verify the dot product of v0 and v1 which is supposed to be 0
#ifdef Debug
        T dot_global = dot(v0_local, v1_local);
        cout << "v"<< iter-1 <<"*v" << iter << " = " << dot_global << endl;
        cout << endl;
        if (std::abs(dot_global) > 1e-5)
            try {
                throw std::runtime_error("Need reorthogonalise: ");
            }
            catch (std::runtime_error& e) {
                std::cerr << "ERROR: " << e.what();
                cout << "v"<< iter-1 <<"*v" << iter << " = " << dot_global << endl;
            }
#endif
    }
    haloUpdate(g_local, v1_local, v1_halo);
    w_local = multGraphVec(g_local, v1_halo);
    alpha_val_global = dot(v1_local, w_local);
    alpha_global.push_back(alpha_val_global);

    transform(g_local, m);
    if (GramSchmidt && g_local.rank() == 0) {
        cout << "Lanczos algorithm WITH GramSchmidt is done." << endl;
    } else if (g_local.rank() == 0) {
        cout << "Lanczos algorithm WITHOUT GramSchmidt is done." << endl;
    }
}
#endif /* end-if Gram-Schmidt */

/*-----------------------------------------------------------------------------
 *  Selective Orthogonalisation
 *-----------------------------------------------------------------------------*/
#ifdef SO_
template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g_local, const int& num_of_eigenvec, bool SO) {
    //VT_TRACER("LANCZOS");
    int local_size = g_local.localSize(), t = 0;
    int global_size = g_local.globalSize();
    int m, scale;
    double tol = 1e-6;

    if (num_of_eigenvec == 1) {
        //SO = false;
        scale = 4 * num_of_eigenvec;
    } else if (num_of_eigenvec == 2) {
        scale = 4 * (num_of_eigenvec - 1);
    } else {
        scale = num_of_eigenvec + 2;
    }
	//cout << "scale = " << scale << endl;
    if (round(log10(global_size)) > 3) {
        scale -= round(log10(std::sqrt(global_size)));
        scale = scale <= 0 ? 1:scale;
    }
    m = scale * std::sqrt(global_size) < global_size ? scale * std::sqrt(global_size):global_size;
	//int size = global_size;
	//cout << "sqrt(" << size << ") = " << std::sqrt(size) << "(" << round(std::sqrt(size)) << "), log10(std::sqrt(" << size << ")) = " << log10(std::sqrt(size)) << "(" << round(log10(std::sqrt(size))) << ")" << endl;
 
	//cout << "scale = " << scale << endl;
	cout << "m = " << m << endl;

    Vector v1_halo(global_size);
    Vector v0_local = init(g_local);
    Vector v1_local = v0_local, w_local, v0_start = v0_local;
    T alpha_val_global = 0.0, beta_val_global = 0.0;

    lanczos_vecs_local[0] = v1_local;
    haloInit(g_local);

    for (int iter = 1; iter < m; iter++) {
        haloUpdate(g_local, v1_local, v1_halo);
        w_local = multGraphVec(g_local, v1_halo);

        alpha_val_global = dot(v1_local, w_local);
        alpha_global.push_back(alpha_val_global);

        for (int i = 0; i < local_size; i++) {
            w_local[i] = w_local[i] - alpha_val_global * v1_local[i] - beta_val_global * v0_local[i];
        }

        beta_val_global = sqrt(dot(w_local, w_local));
        beta_global.push_back(beta_val_global);

        if (std::abs(beta_val_global) < 1e-5) {
            try {
                throw std::runtime_error("Value of beta is close to 0: ");
            }
            catch (std::runtime_error& e) {
                std::cerr << "ERROR: " << e.what();
                cout << "beta[" << iter-1 << "]: " << beta_val_global << endl;
            }
        }

        for (int i = 0; i < local_size; i++) {
            v1_local[i] = w_local[i]/beta_val_global;
        }

        if (std::abs(dot(v0_start, v1_local)) >= tol) {
            gramSchmidt(iter, v1_local);
            t++;
        }
        lanczos_vecs_local[iter] = v1_local;
        v0_local = lanczos_vecs_local[iter-1];
    }
    haloUpdate(g_local, v1_local, v1_halo);
    w_local = multGraphVec(g_local, v1_halo);
    alpha_val_global = dot(v1_local, w_local);
    alpha_global.push_back(alpha_val_global);

    transform(g_local, m);
    if (SO && g_local.rank() == 0) {
        //cout << "m = " << m << endl;
        //cout << "t = " << t << endl;
        cout << "Lanczos algorithm WITH Selective Orthogonalisation is done." << endl;
    } else if (g_local.rank() == 0) {
        //cout << "m = " << m << endl;
        //cout << "t = " << t << endl;
        cout << "Lanczos algorithm WITHOUT Selective Orthogonalisation is done." << endl;
    }
}
#endif /* end-if SO */

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  haloInit
 *  Description:  Find out which rank and the corresponding data need to receive
 * =====================================================================================
 */

template<typename Vector, typename T>
void Lanczos<Vector, T>::haloInit(const Graph& g) {
    // Find out which rank and the corresponding data need to receive
    int size = g.localSize();
    for (int vertex = 0; vertex < size; vertex++) {
        auto it = g.find(g.globalIndex(vertex));
        if (!it->second.empty()) {
            for (const int& neighbour:it->second) {
                int rank = neighbour / size;
                if (rank != g.rank()) {
                    auto it = halo_recv.find(rank);
                    if (it != halo_recv.end()) {
                        it->second.insert(neighbour);
                    } else {
                        std::unordered_set<int> halo_neighbours;
                        halo_neighbours.insert(neighbour);
                        halo_recv.insert(make_pair(rank, halo_neighbours));
                    }
                }
                if (rank != g.rank()) {
                    auto it = halo_send.find(rank);
                    if (it != halo_send.end()) {
                        it->second.insert(g.globalIndex(vertex));
                    } else {
                        std::unordered_set<int> halo_neighbours;
                        halo_neighbours.insert(g.globalIndex(vertex));
                        halo_send.insert(make_pair(rank, halo_neighbours));
                    }
                }
            }
        }
    }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  haloUpdate
 *  Description:  Refresh the halo elements each iteration for Graph * Lanczos_Vec
 * =====================================================================================
 */

template<typename Vector, typename T>
void Lanczos<Vector, T>::haloUpdate(const Graph& g, Vector& v_local, Vector& v_halo) {
    std::vector<mpi::request> reqs;
    std::unordered_map<int, std::unordered_map<int, T>> buf_send; // <global_index, value>;
    std::unordered_map<int, std::unordered_map<int, T>> buf_recv;
    for (int rank = 0; rank < world.size(); rank++) {
        if (rank != g.rank()) {
            auto it = halo_send.find(rank);
            if (it != halo_send.end()) {
                std::unordered_map<int, T> buf_temp;
                for (const int& halo_neighbour:it->second) {
                    buf_temp.insert({halo_neighbour, v_local[g.localIndex(halo_neighbour)]}); // global index, locale value
                }
                buf_send[rank] = buf_temp;
                reqs.push_back(world.isend(rank, 0, buf_send[rank])); //(dest, tag, value to send)
            }
        }
    }
    for (int rank = 0; rank < world.size(); rank++) {
        if (rank != g.rank()) {
            auto it = halo_recv.find(rank);
            if (it != halo_recv.end()) {
                reqs.push_back(world.irecv(rank, 0, buf_recv[rank])); //(src, tag, store to value)
            }
        } else {
            for (int j = 0; j < g.localSize(); j++) {
                v_halo[g.globalIndex(j)] = v_local[j];
            }
        }
    }
    mpi::wait_all(reqs.begin(), reqs.end());
    // Unpack the buffer to fill in to v_halo
    for (int rank = 0; rank < world.size(); rank++) {
        if (rank != g.rank()) {
            std::unordered_map<int, T> buf_temp;
            buf_temp = buf_recv[rank];
            auto it = halo_recv.find(rank);
            if (it != halo_recv.end()) {
                for (const int& halo_neighbour:it->second) {
                    v_halo[halo_neighbour] = buf_temp[halo_neighbour]; //(src, tag, store to value)
                }
            }
        }
    }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  multGraphVec
 *  Description:  The first component of Lanczos iteration fomular, Laplacian matrix * vector
 * =====================================================================================
 */

template<typename Vector, typename T>
Vector Lanczos<Vector, T>::multGraphVec(const Graph& g, const Vector& vec) {
    Vector prod;
    if (g.globalSize() != (int)vec.size())
        throw std::length_error("Lanczos - multGraphVec: The sizes don't match.");

    // Calcualte a partial result in each process, the index starts from zero in vector "prod"
    int size = g.localSize();
    for (int vertex = 0; vertex < size; vertex++) {
        auto it = g.find(g.globalIndex(vertex));
        T temp = 0.0;
        if (!it->second.empty()) {
            for (const int& neighbour:it->second) {
                temp += vec[neighbour];
            }
        }
        prod.push_back(it->second.size() * vec[g.globalIndex(vertex)] - temp);
    }
    return prod;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  gramSchmidt
 *  Description:  GramSchmidt
 * =====================================================================================
 */

template<typename Vector, typename T>
inline void Lanczos<Vector, T>::gramSchmidt(const int& k, Vector& v) {
    //VT_TRACER("GramSchmidt");
    int local_size = lanczos_vecs_local[0].size();
    for (int i = 0; i < k; i++) {
        T dot_global = dot(lanczos_vecs_local[i], v);
        for (int j = 0; j < local_size; j++) {
            v[j] -= dot_global * lanczos_vecs_local[i][j];
        }
    }
    T norm_global = std::sqrt(dot(v, v));
    normalise(v, norm_global);
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  utilities
 *  Description:  Vector operations
 * =====================================================================================
 */

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::dot(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size())
        throw std::length_error("Lanczos - dot: The vector sizes don't match.");

    int local_size = v1.size();
    T dot_local = 0.0, dot_global;
    for (int i = 0; i < local_size; i++) {
        dot_local += v1[i] * v2[i];
    }
    mpi::all_reduce(world, dot_local, dot_global, std::plus<T>());

    return dot_global;
}

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::dot_local(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size())
        throw std::length_error("Lanczos - dot_local: The vector sizes don't match.");

    int local_size = v1.size();
    T dot_local = 0.0;
    for (int i = 0; i < local_size; i++) {
        dot_local += v1[i] * v2[i];	// !!!Index for v1 should be global index.
    }

    return dot_local;
}

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::norm(const Vector& vec) {
    T norm_local = 0.0;
    for (const auto& x:vec) {
        norm_local += x * x;
    }
    return sqrt(norm_local);
}

template<typename Vector, typename T>
inline void Lanczos<Vector, T>::normalise(Vector& vec, const T& norm_global) {
    for (auto x:vec) {
        x /= norm_global;
    }
}

template<typename Vector, typename T>
Vector Lanczos<Vector, T>::init(const Graph& g) {
    int local_size = g.size();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> gen(0.0,1.0);

    Vector vec(local_size);
    for (auto& x:vec) {
        x = gen(generator);
        //x = i;
    }
    T norm_local = norm(vec);
    for (auto& x:vec) {
        x /= norm_local;
        x /= sqrt(g.globalSize()/local_size); // vec[i]/=sqrt(procs), to make sure the global vector is normalised
    }
    return vec;
}

template<typename Vector, typename T>
void Lanczos<Vector, T>::transform(const Graph& g, const int& m) {

    std::vector<std::unordered_map<int, Vector>> gather;
    all_gather(world, lanczos_vecs_local, gather);

    Vector vinitial(g.globalSize(), 0);
    for(int i = 0; i < m; i++)	lanczos_vecs_global[i] = vinitial;

    for (int row = 0; row < m; row ++) {
        for (int i = 0; i < world.size(); i++) {
            for (int j = 0; j < g.localSize(); j++) {
                lanczos_vecs_global[row][i * g.localSize() + j] = gather[i][row][j];
            }
        }
    }
}

template<typename Vector, typename T>
void Lanczos<Vector, T>::print_tri_mat() {
    int size = alpha_global.size();
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            if (row == col)
                cout << alpha_global[row] << "\t";
            else if (col - row == 1)
                cout << beta_global[row] << "\t";
            else if (row - col == 1)
                cout << beta_global[col] << "\t";
            else
                cout << "0" << "\t";
        }
        cout << endl;
    }
}
#endif
