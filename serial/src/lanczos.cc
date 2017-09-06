/**
 * @file lanczos.cc
 * @brief Lanczos algorithm
 * @author Ken Hu, xnchnhu@gmail.com
 */

#ifndef LANCZOS_CC_
#define LANCZOS_CC_

#include "lanczos.h"
#include <cmath>
#include <exception>
#include <iostream>
#include <utility>

#ifdef VT_
#include "vt_user.h"
#endif

using std::cout;
using std::endl;

/**
 * @brief Lanczos algorithm with selective orthogonalisation
 * @param FILL-ME-IN
 * @return FILL-ME-IN
 */

/*-----------------------------------------------------------------------------
 * 	input:
 *		G contains all the elements of the original graph
 *
 *	output:
 *		alpha[0..n-1] returns all the diagonal elements of the
 *tridiagonal
 *matrix
 *		beta[0..n-2] returns all the subdiagonal elements of the
 *tridiagonal matrix
 *		lanczos_vecs[0..n-1][0..n-1] returns n vectors, each row
 *represents a lanczos vector
 *-----------------------------------------------------------------------------*/

#include <cmath>
template <typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g, const int& num_of_eigenvec, bool SO)
{
#ifdef VT_
    VT_TRACER("LANCZOS_SO");
#endif
    const int size = g.size();
    int t = 0;
    int m = getIteration(num_of_eigenvec, size);

    Vector v0 = init(size);
    Vector v1 = v0, w, vstart = v0;

    T beta_val = 0.0, tol = 1e-6;
    alpha.resize(m);
    beta.resize(m - 1);
    lanczos_vecs.resize(m);
    lanczos_vecs[0] = v0;

    for (int iter = 1; iter < m; iter++) {
        w = multGraphVec(g, v1);
        alpha[iter - 1] = dot(v1, w);
        for (int i = 0; i < size; i++) {
            w[i] = w[i] - alpha[iter - 1] * v1[i] - beta_val * v0[i];
        }

        beta_val = norm(w);
        beta[iter - 1] = beta_val;
        /*
        if (std::abs(beta[iter - 1]) < 1e-5) {
            try {
                throw std::runtime_error("Value of beta is close to 0: ");
            }
            catch (std::runtime_error& e) {
                std::cerr << "ERROR: " << e.what();
                cout << "beta[" << iter - 1 << "]: " << beta[iter - 1] << endl;
            }
        }
        */
        for (int index = 0; index < size; index++) {
            v1[index] = w[index] / beta[iter - 1];
        }
        if (SO) {
            if (std::abs(dot(vstart, v1)) >= tol) {
                gramSchmidt(iter, v1);
                t++;
            }
        }
        lanczos_vecs[iter] = v1;
        v0 = lanczos_vecs[iter - 1];
    }
    w = multGraphVec(g, v1);
    alpha[m - 1] = dot(v1, w);
    if (SO) {
        cout << "Lanczos algorithm WITH Selective Orthogonalisation is done."
             << endl;
    } else {
        cout << "Lanczos algorithm WITHOUT Selective Orthogonalisation is done."
             << endl;
    }
    cout << "number of iterations = " << m
         << ", number of Orthogonalisation = " << t << endl;
}

/**
 * @brief Calculate iterations for Lanczos algorithm
 * @param FILL-ME-IN
 * @return FILL-ME-IN
 */

template <typename Vector, typename T>
const int Lanczos<Vector, T>::getIteration(const int& num_of_eigenvec,
                                           const int& size)
{
    int scale, m;
    if (num_of_eigenvec == 1) {
        scale = 4 * num_of_eigenvec;
    } else if (num_of_eigenvec == 2) {
        scale = 4 * (num_of_eigenvec - 1);
    } else {
        scale = num_of_eigenvec + 2;
    }
    if (round(log10(size)) > 3) {
        scale -= round(log10(std::sqrt(size)));
        scale = scale <= 0 ? 1 : scale;
    }
    m = scale * std::sqrt(size) < size ? scale * std::sqrt(size) : size;
    return m;
}

/**
 * @brief The first component of Lanczos iteration fomular, Laplacian matrix *
 * vector
 * @param FILL-ME-IN
 * @return FILL-ME-IN
 */

template <typename Vector, typename T>
Vector Lanczos<Vector, T>::multGraphVec(const Graph& g, const Vector& vec)
{
    Vector prod(g.size());
    for (auto it = g.cbegin(); it != g.cend(); it++) {
        T temp = 0.0;
        if (!it->second.empty()) {
            for (const int& neighbour : it->second) {
                temp += vec[neighbour];
            }
        }
        prod[it->first] = it->second.size() * vec[it->first] - temp;
    }
    return prod;
}

/**
 * @brief Reorthogonalisation
 * @param FILL-ME-IN
 * @return FILL-ME-IN
 */

template <typename Vector, typename T>
inline void Lanczos<Vector, T>::gramSchmidt(const int& k, Vector& v)
{
#ifdef VT_
    VT_TRACER("GramSchmidt");
#endif
    int size = v.size();
    for (int i = 0; i < k; i++) {
        T reorthog_dot_product = dot(lanczos_vecs[i], v);
        // cout << "iter " << i << " gramSchmidt global dot product " <<
        // reorthog_dot_product << endl;
        for (int j = 0; j < size; j++) {
            v[j] -= reorthog_dot_product * lanczos_vecs[i][j];
        }
    }
    normalise(v);
}

/**
 * @brief Vector operations
 * @param FILL-ME-IN
 * @return FILL-ME-IN
 */

template <typename Vector, typename T>
inline T Lanczos<Vector, T>::dot(const Vector& v1, const Vector& v2)
{
    if (v1.size() != v2.size())
        throw std::length_error("Lanczos - dot: The vector sizes don't match.");
    int size = v1.size();
    T dotprod = 0.0;
    for (int index = 0; index < size; index++) {
        dotprod += v1[index] * v2[index];
    }
    return dotprod;
}

template <typename Vector, typename T>
inline T Lanczos<Vector, T>::norm(const Vector& vec)
{
    T normret = 0.0;
    for (const auto& x : vec) {
        normret += x * x;
    }
    return sqrt(normret);
}

template <typename Vector, typename T>
inline T Lanczos<Vector, T>::l2norm(const Vector& alpha, const Vector& beta)
{
    T normret = 0.0;
    T col_sum = 0.0;

    int size = alpha.size();
    for (int col = 0; col < size; col++) {
        col_sum = 0.0;
        for (int row = 0; row < size; row++) {
            if (row == col)
                col_sum = alpha[row] * alpha[row];
            else if (col - row == 1)
                col_sum += beta[row] * beta[row];
            else if (row - col == 1)
                col_sum += beta[col] * beta[col];
        }
        normret += std::sqrt(col_sum);
    }
    cout << "norm = " << normret << endl;
    return normret;
}

template <typename Vector, typename T>
inline Vector& Lanczos<Vector, T>::normalise(Vector& vec)
{
    T normal = norm(vec);
    for (auto& x : vec) {
        x /= normal;
    }
    return vec;
}

template <typename Vector, typename T>
Vector Lanczos<Vector, T>::init(const int& size)
{
    Vector vec(size);
    for (auto& x : vec) {
        x = drand48();
    }
    T normalise = norm(vec);
    for (auto& x : vec) {
        x /= normalise;
    }
    return vec;
}

template <typename Vector, typename T>
void Lanczos<Vector, T>::print_tri_mat()
{
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
                cout << "0"
                     << "\t";
        }
        cout << endl;
    }
}
#endif
