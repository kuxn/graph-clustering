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
#include <set>

#include <boost/serialization/serialization.hpp>
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
    int local_size = g_local.size();
    int m = getIteration(num_of_eigenvec, global_size);
    Vector v1_halo(g_local.globalSize());

    Vector v0_local = init(g_local);
    Vector v1_local = v0_local, w_local;
    T alpha_val_global = 0.0, beta_val_global = 0.0;

    lanczos_vecs[0] = v1_local;
    haloInit(g_local);

    for (int iter = 1; iter < m; iter++) {
        haloUpdate(g_local, v1_local, v1_halo);
        w_local = multGraphVec(g_local, v1_halo);
        alpha_val_global = dot(v1_local, w_local);
        alpha.push_back(alpha_val_global);

        for (int i = 0; i < local_size; i++) {
            w_local[i] = w_local[i] - alpha_val_global * v1_local[i] - beta_val_global * v0_local[i];
        }

        beta_val_global = sqrt(dot(w_local, w_local));
        beta.push_back(beta_val_global);

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
        if (GramSchmidt) {
            gramSchmidt(iter, v1_local);
        }
        lanczos_vecs[iter] = v1_local;
        v0_local = lanczos_vecs[iter-1];

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
    alpha.push_back(alpha_val_global);

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
    int local_size = g_local.size();
    int global_size = g_local.globalSize();
    int m, t = 0;
    double tol = 1e-6;
    m = getIteration(num_of_eigenvec, global_size);

    Vector v1_halo(global_size);
    Vector v0_local = init(g_local);

    Vector v1_local = v0_local, w_local, v0_start = v0_local;
    T alpha_val_global = 0.0, beta_val_global = 0.0;

    lanczos_vecs.resize(m);
    lanczos_vecs[0] = v1_local;
    haloInit(g_local);

    for (int iter = 1; iter < m; iter++) {
        haloUpdate(g_local, v1_local, v1_halo);
        w_local = multGraphVec(g_local, v1_halo);
        alpha_val_global = dot(v1_local, w_local);
        alpha.push_back(alpha_val_global);

        for (int i = 0; i < local_size; i++) {
            w_local[i] = w_local[i] - alpha_val_global * v1_local[i] - beta_val_global * v0_local[i];
        }
        beta_val_global = sqrt(dot(w_local, w_local));
        beta.push_back(beta_val_global);

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
        if (SO && std::abs(dot(v0_start, v1_local)) >= tol) {
            gramSchmidt(iter, v1_local);
            t++;
        }

        lanczos_vecs[iter] = v1_local;
        v0_local = lanczos_vecs[iter-1];
    }
    haloUpdate(g_local, v1_local, v1_halo);
    w_local = multGraphVec(g_local, v1_halo);
    alpha_val_global = dot(v1_local, w_local);
    alpha.push_back(alpha_val_global);

    if (g_local.rank() == 0) {
        cout << "m = " << m << ", t = " << t << endl;
    }
    if (SO && g_local.rank() == 0) {
        cout << "Lanczos algorithm WITH Selective Orthogonalisation is done." << endl;
    } else if (g_local.rank() == 0) {
        cout << "Lanczos algorithm WITHOUT Selective Orthogonalisation is done." << endl;
    }
}
#endif /* end-if SO */

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  getIteration
 *  Description:  Calculate iterations for Lanczos algorithm
 * =====================================================================================
 */

template<typename Vector, typename T>
const int Lanczos<Vector, T>::getIteration(const int& num_of_eigenvec, const int & global_size) {
    int scale, m;
    if (num_of_eigenvec == 1) {
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
    m = global_size;
    return m;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  haloInit
 *  Description:  Find out which rank and the corresponding data need to receive
 * =====================================================================================
 */

template<typename Vector, typename T>
void Lanczos<Vector, T>::haloInit(const Graph& g) {
    // Find out which rank and the corresponding data need to receive
    std::unordered_map<int, std::set<int>> halo_recv_temp; // <rank, halo_neighbours to receive>
    std::unordered_map<int, std::set<int>> halo_send_temp; // <rank, halo_neighbours to send>
    for (auto iter = g.cbegin(); iter != g.cend(); ++iter) {
        if (!iter->second.empty()) {
            for (const int& neighbour:iter->second) {
                int rank = g.global_rank_map[neighbour];
                if (rank != g.rank()) {
                    auto it = halo_recv_temp.find(rank);
                    if (it != halo_recv_temp.end()) {
                        it->second.insert(neighbour);
                    } else {
                        std::set<int> halo_neighbours;
                        halo_neighbours.insert(neighbour);
                        halo_recv_temp.insert({rank, halo_neighbours});
                    }
                }
                if (rank != g.rank()) {
                    auto it = halo_send_temp.find(rank);
                    if (it != halo_send_temp.end()) {
                        it->second.insert(iter->first);
                        //it->second.insert(g.globalIndex(vertex));
                    } else {
                        std::set<int> halo_neighbours;
                        halo_neighbours.insert(iter->first);
                        //halo_neighbours.insert(g.globalIndex(vertex));
                        halo_send_temp.insert({rank, halo_neighbours});
                    }
                }
            }
        }
    }
    // Convert <int, set> to <int, vector> for looking up in halo_update, cheaper than iterating a set.
    for (auto& it:halo_send_temp) {
        int rank = it.first;
        std::vector<int> vector_send;
        for (auto& x:it.second) {
            vector_send.push_back(x);
        }
        halo_send.insert({rank, vector_send});
    }
    for (auto& it:halo_recv_temp) {
        int rank = it.first;
        std::vector<int> vector_recv;
        for (auto& x:it.second) {
            vector_recv.push_back(x);
        }
        halo_recv.insert({rank, vector_recv});
    }
    //if (g.rank() == 2) {
    //	cout << "in rank " << g.rank() << " ";
    //	for (auto& it:halo_send) {
    //		cout << "send to rank " << it.first << ": " << endl;
    //		for (unsigned int i = 0; i < it.second.size(); ++i) {
    //			cout << "vector_send["<< i << "] " << it.second[i] << endl;
    //		}
    //	}
    //	cout << "in rank " << g.rank() << " ";
    //	for (auto& it:halo_recv) {
    //		cout << "recv from rank " << it.first << ": " << endl;
    //		for (unsigned int i = 0; i < it.second.size(); ++i) {
    //			cout << "vector_recv["<< i << "] " << it.second[i] << endl;
    //		}
    //	}
    //}
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
    std::unordered_map<int, std::vector<T>> buf_send(world.size()); // <global_index, value>;
    std::unordered_map<int, std::vector<T>> buf_recv(world.size());

    for (int rank = 0; rank < world.size(); rank++) {
        if (rank != g.rank()) {
            auto it = halo_send.find(rank);
            if (it != halo_send.end()) {
                std::vector<T> buf_temp;
                for (const int& halo_neighbour:it->second) {
                    buf_temp.push_back(v_local[g.localIndex(halo_neighbour)]);
                }
                buf_send.insert({rank, buf_temp});
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
            for (int j = 0; j < g.size(); j++) {
                v_halo[g.globalIndex(j)] = v_local[j];
            }
        }
    }
    mpi::wait_all(reqs.begin(), reqs.end());
    // Unpack the buffer to fill in to v_halo
    for (int rank = 0; rank < world.size(); rank++) {
        if (rank != g.rank()) {
            std::vector<T> buf_temp;
            buf_temp = buf_recv[rank];
            auto it = halo_recv.find(rank);
            if (it != halo_recv.end()) {
                int i = 0;
                for (const int& halo_neighbour:it->second) {
                    v_halo[halo_neighbour] = buf_temp[i]; //(src, tag, store to value)
                    i++;
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
    Vector prod(g.size());
    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        T temp = 0.0;
        if (!it->second.empty()) {
            for (const int& neighbour:it->second) {
                temp += vec[neighbour];
            }
        }
        prod[g.localIndex(it->first)] = it->second.size() * vec[it->first] - temp;
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
    int local_size = lanczos_vecs[0].size();
    for (int i = 0; i < k; i++) {
        T dot_global = dot(lanczos_vecs[i], v);
        for (int j = 0; j < local_size; j++) {
            v[j] -= dot_global * lanczos_vecs[i][j];
        }
    }
    // Normalise
    T norm_global = std::sqrt(dot(v, v));
    for (auto& x:v) {
        x /= norm_global;
    }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  utilities
 *  Description:  Vector operations
 * =====================================================================================
 */

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::dot(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size()) {
        throw std::length_error("Lanczos - dot: The vector sizes don't match.");
    }
    int local_size = v1.size();
    T dot_local = 0.0, dot_global;
    for (int i = 0; i < local_size; i++) {
        dot_local += v1[i] * v2[i];
    }
    mpi::all_reduce(world, dot_local, dot_global, std::plus<T>());

    return dot_global;
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
Vector Lanczos<Vector, T>::init(const Graph& g) {
    int local_size = g.size();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> gen(0.0,1.0);

    Vector vec(local_size);
    for (auto& x:vec) {
        x = gen(generator);
    }
    T norm_local = norm(vec);
    for (auto& x:vec) {
        x /= norm_local;
        x /= sqrt(world.size()); // vec[i]/=sqrt(procs), to make sure the global vector is normalised
    }
    return vec;
}

template<typename Vector, typename T>
void Lanczos<Vector, T>::print_tri_mat() {
    int size = alpha.size();
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            if (row == col)
                cout << alpha[row] << "\t";
            else if (col - row == 1)
                cout << beta[row] << "\t";
            else if (row - col == 1)
                cout << beta[col] << "\t";
            else
                cout << "0" << "\t";
        }
        cout << endl;
    }
}
#endif
