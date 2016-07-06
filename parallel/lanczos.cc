/*
 * =====================================================================================
 *
 *       Filename:  lanczos.cpp
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

#include "lanczos.h"
#include "vt_user.h"

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
 *  Modified Lanczos algorithm with reorthogonalisation by Gram–Schmidt
 *-----------------------------------------------------------------------------*/

template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(mpi::communicator& world, const Graph& g_local, bool reorthogonalisation) {

    VT_TRACER("LANCZOS");
    int local_size = g_local.localSize();
    Vector v0_local(local_size);
    v0_local = init(v0_local, g_local.globalSize());

#ifdef Debug
    cout << "rank " << g_local.rank() << ": ";
    for (int i = 0; i < local_size; i++) {
        cout << "v0_local[" << i << "] = " << v0_local[i] << " ";
    }
    cout << endl;
    cout << "global dot of these vectors = " << dot(world, v0_local, v0_local) << endl;
#endif

    Vector t_local = v0_local, v1_local = v0_local, w_local;

    T alpha_val_global = 0.0, beta_val_global = 0.0;

    std::vector<Vector> v0_global_gather, v1_global_gather;
    Vector v0_global, v1_global;
    mpi::all_gather(world, v0_local, v0_global_gather); // Gather all local vectors to global
    v0_global = transform(v0_global_gather, local_size);
    v1_global = v0_global;

#ifdef Debug
    cout << "v1_global:";
    for (auto x:v1_global)
        cout << x << " ";
#endif

    lanczos_vecs_global[0] = v0_global;

    for (int iter = 1; iter < g_local.globalSize(); iter++) {
        w_local = multGraphVec(g_local, v1_global);

#ifdef Debug
        cout << "w_local, rank " << g_local.rank() << ": ";
        for (auto x:w_local)
            cout << x << " ";
        cout << endl;
#endif

        T alpha_val_global = dot(world, v1_local, w_local);
        //cout << "dot(v1, v1) = " << dot(v1, v1) << endl;
        alpha_global.push_back(alpha_val_global);

#ifdef Debug
        cout << "rank " << g_local.rank() << ": ";
        cout << "alpha[" << iter - 1 << "] =  " << alpha_val_global << endl;
#endif

        for (int i = 0; i < local_size; i++) {
            t_local[i] = w_local[i] - alpha_val_global * v1_local[i] - beta_val_global * v0_local[i];
        }

        beta_val_global = sqrt(dot(world, t_local, t_local)); 
        beta_global.push_back(beta_val_global);	

        if (std::abs(beta_val_global) < 1e-5) 
            try { throw std::runtime_error("Value of beta is close to 0: "); }
        catch (std::runtime_error& e) { 
            std::cerr << "ERROR: " << e.what(); 
            cout << "beta[" << iter-1 << "]: " << beta_val_global << endl;
        }

        if (!reorthogonalisation) {
            v0_local = v1_local;
        }
        for (int i = 0; i < local_size; i++) {
            v1_local[i] = t_local[i]/beta_val_global;
        }

        mpi::all_gather(world, v1_local, v1_global_gather); // Gather all local vectors to global
        v1_global = transform(v1_global_gather, local_size);
        lanczos_vecs_global[iter] = v1_global;

        if (reorthogonalisation) {
            gramSchmidt(iter, g_local.globalSize());
            v0_global = lanczos_vecs_global[iter-1];
            v1_global = lanczos_vecs_global[iter];
            for (int i = 0; i < local_size; i++) {
                v0_local[i] = v0_global[g_local.globalIndex(i)];
                v1_local[i] = v1_global[g_local.globalIndex(i)];
            }
        }

        //Verify the dot product of v0 and v1 which is supposed to be 0
        T dot_global = dot(world, v0_local, v1_local);
#ifdef Debug
        cout << "v"<< iter-1 <<"*v" << iter << " = " << dot_global << endl;
        cout << endl;
#endif
        if (std::abs(dot_global) > 1e-5) 
            try { throw std::runtime_error("Need reorthogonalise: "); }
        catch (std::runtime_error& e) { 
            std::cerr << "ERROR: " << e.what(); 
            cout << "v"<< iter-1 <<"*v" << iter << " = " << dot_global << endl;
        }
    }
    w_local = multGraphVec(g_local, v1_global);
    alpha_val_global = dot(world, v1_local, w_local);
    alpha_global.push_back(alpha_val_global);

    if (reorthogonalisation) {
        cout << "Lanczos algorithm WITH reorthogonalisation is done." << endl;
    } else {
        cout << "Lanczos algorithm WITHOUT reorthogonalisation is done." << endl; 
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
    //if (g.size() != (int)vec.size())
    //throw std::length_error("The sizes don't match.");

    // Calcualte a partial result in each process, the index starts from zero in vector "prod"

    int size = g.localSize();
    for (int vertex = 0; vertex < size; vertex++) {
        auto it = g.find(g.globalIndex(vertex));
        T temp = 0.0;
        if (!it->second.empty())
            for (const int& neighbour:it->second) {
                temp += vec[neighbour];
                //if (g.rank() == 3)
                //cout << "rank " << g.rank() << " vec[" << neighbour << "] = " << vec[neighbour] << endl;
            }

        //if (g.rank() == 3)
        //cout << "push_back " << it->second.size() << " * " << vec[g.globalIndex(vertex)] << endl;
        prod.push_back(it->second.size() * vec[g.globalIndex(vertex)] - temp);
    }
    return prod;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  gramSchmidt
 *  Description:  Reorthogonalisation
 * =====================================================================================
 */

template<typename Vector, typename T>
inline void Lanczos<Vector, T>::gramSchmidt(int& iter, const int& size) {
    VT_TRACER("GramSchmidt");
    for (int k = 1; k <= iter; k++) {
        //cout << "i - norm of lanczos_vecs_global["<<k<<"] = " << norm(lanczos_vecs_global[k]) << endl;
        for (int i = 0; i < k; i++) {
            T reorthog_dot_product = dot_local(lanczos_vecs_global[i], lanczos_vecs_global[k]);
            for (int j = 0; j < size; j++) {
                lanczos_vecs_global[k][j] -= reorthog_dot_product * lanczos_vecs_global[i][j];
            }
        }
        //T normalise = norm(lanczos_vecs_global[k]);
        //for (int j = 0; j < size; j++) lanczos_vecs_global[k][j] /= normalise;
        normalise(lanczos_vecs_global[k]);
        //cout << "norm of lanczos_vecs_global["<<iter<<"] = " << norm(lanczos_vecs_global[iter]) << endl;
    }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  utilities
 *  Description:  Vector operations
 * =====================================================================================
 */

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::dot(mpi::communicator& world, const Vector& v1, const Vector& v2) {
    //if (v1.size() != v2.size())	
    //throw std::length_error("The vector sizes don't match.");

    int local_size = v1.size();
    T dot_local = 0.0, dot_global;
    for (int i = 0; i < local_size; i++) {
        dot_local += v1[i] * v2[i];	// !!!Index for v1 should be global index.
    }

    mpi::all_reduce(world, dot_local, dot_global, std::plus<T>());

    return dot_global;
}

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::dot_local(const Vector& v1, const Vector& v2) {
    //if (v1.size() != v2.size())	
    //throw std::length_error("The vector sizes don't match.");

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
    for (const T& value:vec) {
        norm_local += value * value;
    }
    return sqrt(norm_local);
}

template<typename Vector, typename T>
inline Vector& Lanczos<Vector, T>::normalise(Vector& vec) {
    int local_size = vec.size();
    T norm_local = norm(vec);
    for (int i = 0; i < local_size; i++) {
        vec[i] /= norm_local;
    }
    return vec;
}

template<typename Vector, typename T>
Vector& Lanczos<Vector, T>::init(Vector& vec, int global_size) {
    int local_size = vec.size();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> gen(0.0,1.0);
    for (int i = 0; i < local_size; i++) {
        vec[i] = gen(generator);
    }
    T norm_local = norm(vec);
    for (int i = 0; i < local_size; i++) {
        vec[i] /= norm_local;
        vec[i] /= sqrt(global_size/local_size); // vec[i]/=sqrt(procs), to make sure the global vector is normalised
    }
    return vec;
}

template<typename Vector, typename T>
Vector Lanczos<Vector, T>::transform(std::vector<Vector>& vec_gather, int local_size) {
    Vector global_vec;

    for (unsigned int i = 0; i < vec_gather.size(); i++) {
        for (int j = 0; j < local_size; j++) {
            global_vec.push_back(vec_gather[i][j]);
        }	
    }
    return global_vec;
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
