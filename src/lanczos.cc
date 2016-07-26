/*
 * =====================================================================================
 *
 *       Filename:  lanczos.cpp
 *
 *    Description:  Lanczos algorithm
 *        Created:  06/11/2016 15:33:49
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef LANCZOS_CPP_
#define LANCZOS_CPP_

#include <iostream>
#include <exception>
#include <utility>
#include <cmath>

#include "lanczos.h"
#include "tqli.h"
#include "vt_user.h"

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
 *		beta[0..n-2] returns all the subdiagonal elements of the tridiagonal matrix
 *		lanczos_vecs[0..n-1][0..n-1] the kth row returns the kth Lanczos vector
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 *  Modified Lanczos algorithm with GramSchmidt by Gram–Schmidt
 *-----------------------------------------------------------------------------*/
#ifdef GS_
template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g, const int& num_of_eigenvec, bool GramSchmidt) {
    VT_TRACER("LANCZOS");
    const int size = g.size();
    //const int m = 6 * std::sqrt(size);
    const int m = size;

    Vector v0 = init(size);
    Vector v1 = v0, w;

    T beta_val = 0.0;
    alpha.resize(m);
    beta.resize(m - 1);
    lanczos_vecs[0] = v0;

    for (int iter = 1; iter < m; iter++) {
        w = multGraphVec(g, v1);
        alpha[iter - 1] = dot(v1, w);
        for (int i = 0; i < size; i++) {
            w[i] = w[i] - alpha[iter - 1] * v1[i] - beta_val * v0[i];
        }

        beta_val = norm(w);
        beta[iter - 1] = beta_val;
        if (std::abs(beta[iter - 1]) < 1e-5) {
            try {
                throw std::runtime_error("Value of beta is close to 0: ");
            }
            catch (std::runtime_error& e) {
                std::cerr << "ERROR: " << e.what();
                cout << "beta[" << iter - 1 << "]: " << beta[iter - 1] << endl;
            }
        }

        for (int index = 0; index < size; index++) {
            v1[index] = w[index]/beta[iter - 1];
        }

        if (GramSchmidt) {
            gramSchmidt(iter, v1);
        }
        lanczos_vecs[iter] = v1;
        v0 = lanczos_vecs[iter - 1];

        //Verify the dot product of v0 and v1 which is supposed to be 0
        T dot_product = dot(v0, v1);
        if (std::abs(dot_product) > 1e-5) {
            try {
                throw std::runtime_error("Need reorthogonalise: ");
            }
            catch (std::runtime_error& e) {
                std::cerr << "ERROR: " << e.what();
                cout << "v"<< iter - 1 <<"*v" << iter << " = " << dot_product << endl;
            }
        }
    }
    w = multGraphVec(g, v1);
    alpha[m - 1] = dot(v1, w);

    if (GramSchmidt) {
        cout << "Lanczos algorithm WITH GramSchmidt is done." << endl;
    } else {
        cout << "Lanczos algorithm WITHOUT GramSchmidt is done." << endl;
    }
}
#endif // endif - GS_
/*-----------------------------------------------------------------------------
 *  Modified Lanczos algorithm with selective orthogonalisation
 *-----------------------------------------------------------------------------*/
#ifdef SO_
#include <cmath>
template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g, const int& num_of_eigenvec, bool SO) {
    VT_TRACER("LANCZOS_SO");
    const int size = g.size();
    //const int size = 1e9;
	int m, scale;
	if (num_of_eigenvec == 1) {
		//SO = false;
		scale = 4 * num_of_eigenvec;
	} else if (num_of_eigenvec == 2){
		scale = 4 * (num_of_eigenvec - 1);
	} else {
		scale = num_of_eigenvec + 2;
	}

	cout << "scale = " << scale << endl;
	cout << "sqrt(" << size << ") = " << std::sqrt(size) << "(" << round(std::sqrt(size)) << "), log10(std::sqrt(" << size << ")) = " << log10(std::sqrt(size)) << "(" << round(log10(std::sqrt(size))) << ")" << endl;

	if (round(log10(size)) > 3) {
		scale -= round(log10(std::sqrt(size)));
		scale = scale <= 0 ? 1:scale;
	}

	cout << "scale = " << scale << endl;
	m = scale * std::sqrt(size) < size ? scale *  std::sqrt(size):size;
	m = size;
	cout << "m = " << m << endl;

    Vector v0 = init(size);
    Vector v1 = v0, w, vstart = v0;

    T beta_val = 0.0, tol = 1e-6;
    alpha.resize(m);
    beta.resize(m - 1);
    lanczos_vecs[0] = v0;

    int t = 0;
    for (int iter = 1; iter < m; iter++) {
        w = multGraphVec(g, v1);
        alpha[iter - 1] = dot(v1, w);
        for (int i = 0; i < size; i++) {
            w[i] = w[i] - alpha[iter - 1] * v1[i] - beta_val * v0[i];
        }

        beta_val = norm(w);
        beta[iter - 1] = beta_val;
        //if (std::abs(beta[iter - 1]) < 1e-5) {
        //    try {
        //        throw std::runtime_error("Value of beta is close to 0: ");
        //    }
        //    catch (std::runtime_error& e) {
        //        std::cerr << "ERROR: " << e.what();
        //        cout << "beta[" << iter - 1 << "]: " << beta[iter - 1] << endl;
        //    }
        //}

        for (int index = 0; index < size; index++) {
            v1[index] = w[index]/beta[iter - 1];
        }

        if (SO) {
            //std::unordered_map<int, Vector> q;
            //Vector d = alpha;
            //Vector e = beta;
            //tqli(d, e, q);
            //cout << "beta_val * std::abs(q[iter][iter - 1] = " << beta_val * std::abs(q[iter][iter - 1]) << endl;
            if (std::abs(dot(vstart, v1)) >= tol) {
            //if (beta_val * std::abs(q[iter][iter - 1]) <= tol) {
                gramSchmidt(iter, v1);
                t++;
            }
        }
        lanczos_vecs[iter] = v1;
        v0 = lanczos_vecs[iter - 1];

        //Verify the dot product of v0 and v1 which is supposed to be 0
        //T dot_product = dot(vstart, v1);
        ////cout << "v0 * v1 = " << dot(v0, v1) << endl;
        //if (std::abs(dot_product) > 1e-5) {
        //    try {
        //        throw std::runtime_error("Need reorthogonalise: ");
        //    }
        //    catch (std::runtime_error& e) {
        //        std::cerr << "ERROR: " << e.what();
        //        cout << "v"<< iter - 1 <<"*v" << iter << " = " << dot_product << endl;
        //    }
        //}
    }
    w = multGraphVec(g, v1);
    alpha[m - 1] = dot(v1, w);

    cout << "t = " << t << endl;
    if (SO) {
        cout << "Lanczos algorithm WITH Selective Orthogonalisation is done." << endl;
    } else {
        cout << "Lanczos algorithm WITHOUT Selective Orthogonalisation is done." << endl;
    }
}
#endif // endif - SO

/*-----------------------------------------------------------------------------
 *  Modified Lanczos algorithm with restart
 *-----------------------------------------------------------------------------*/
#ifdef RS_
template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g, const int& num_of_eigenvec, bool RS) {
    VT_TRACER("LANCZOS_RS");
    const int size = g.size();
    const int m = 2 * std::sqrt(size);
    //int m = size;
    Vector v0 = init(size);

    const int MAXIT = 5;
    for (int k = 1; k < MAXIT; k++) {

        Vector v1 = v0, w;
        T beta_val = 0.0, tol = 1e-6;
        alpha.resize(m);
        beta.resize(m - 1);
        lanczos_vecs[0] = v0;

        for (int iter = k; iter < m; iter++) {
            w = multGraphVec(g, v1);
            alpha[iter - 1] = dot(v1, w);
            for (int i = 0; i < size; i++) {
                w[i] = w[i] - alpha[iter - 1] * v1[i] - beta_val * v0[i];
            }
            beta_val = norm(w);
            beta[iter - 1] = beta_val;
            if (std::abs(beta[iter - 1]) < 1e-5) {
                try {
                    throw std::runtime_error("Value of beta is close to 0: ");
                }
                catch (std::runtime_error& e) {
                    std::cerr << "ERROR: " << e.what();
                    cout << "beta[" << iter - 1 << "]: " << beta[iter - 1] << endl;
                }
            }
            for (int index = 0; index < size; index++) {
                v1[index] = w[index]/beta[iter - 1];
            }
            lanczos_vecs[iter] = v1;
            v0 = lanczos_vecs[iter - 1];

            //Verify the dot product of v0 and v1 which is supposed to be 0
            T dot_product = dot(v0, v1);
            if (std::abs(dot_product) > 1e-5) {
                try {
                    throw std::runtime_error("Need reorthogonalise: ");
                }
                catch (std::runtime_error& e) {
                    std::cerr << "ERROR: " << e.what();
                    cout << "v"<< iter - 1 <<"*v" << iter << " = " << dot_product << endl;
                }
            }
        }
        w = multGraphVec(g, v1);
        alpha[m - 1] = dot(v1, w);
        std::unordered_map<int, Vector> q;
        Vector d = alpha;
        Vector e = beta;
        tqli(d, e, q);
        // Calculate eigenvectors
        cout << endl;
        std::unordered_map<int, Vector> ritz_vector;
        Vector vinit(size, 0);
        for(int i = 0; i < m; i++)	ritz_vector[i] = vinit;
        for (int row = 0; row < m; row++) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < m; j++) {
                    ritz_vector[row][i] += lanczos_vecs[j][i] * q[j][row];
                }
            }
        }

        cout << "Ritz values: " << endl;
        for (auto& x:d) {
            cout << x << " ";
        }
        cout << endl;
        cout << "Ritz vectors: " << endl;
        int row_size = ritz_vector.size();
        for (int row = 0; row < row_size; row++) {
            int col_size = ritz_vector[row].size();
            for (int col = 0; col < col_size; col++) {
                cout << ritz_vector[row][col] << " ";
            }
            cout << endl;
        }
        // new v0 for restart
        T temp = 0.0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < m; j++) {
                temp += ritz_vector[j][i];
            }
            v0[i] = temp;
        }
        normalise(v0);
        cout << "norm of v0 = " << norm(v0) << endl;

        //if (beta_val * std::abs(q[iter][iter - 1]) <= tol) {
        //	gramSchmidt(iter, v1);
        //	t++;
        //}
    }

    //cout << "t = " << t << endl;
    if (RS) {
        cout << "Lanczos algorithm WITH Restart is done." << endl;
    } else {
        cout << "Lanczos algorithm WITHOUT Restart is done." << endl;
    }
}
#endif // endif - RS

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  multGraphVec
 *  Description:  The first component of Lanczos iteration fomular, Laplacian matrix * vector
 * =====================================================================================
 */

template<typename Vector, typename T>
Vector Lanczos<Vector, T>::multGraphVec(const Graph& g, const Vector& vec) {
    Vector prod;
    if (g.size() != (int)vec.size())
        throw std::length_error("Lanczos - multGraphVec: The sizes don't match.");

    int numofvertex = vec.size();
    for (int vertex = 0; vertex < numofvertex; vertex++) {
        auto it = g.find(vertex);
        T temp = 0.0;
        if (!it->second.empty()) {
            for (const int& neighbour:it->second) {
                temp += vec[neighbour];
            }
        }
        prod.push_back(it->second.size() * vec[vertex] - temp);
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
inline void Lanczos<Vector, T>::gramSchmidt(const int& k,  Vector& v) {
    VT_TRACER("GramSchmidt");
    int size = v.size();
    for (int i = 0; i < k; i++) {
        T reorthog_dot_product = dot(lanczos_vecs[i], v);
        //cout << "iter " << i << " gramSchmidt global dot product " << reorthog_dot_product << endl;
        for (int j = 0; j < size; j++) {
            v[j] -= reorthog_dot_product * lanczos_vecs[i][j];
        }
    }
    normalise(v);
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
    int size = v1.size();
    T dotprod = 0.0;
    for (int index = 0; index < size; index++) {
        dotprod += v1[index] * v2[index];
    }
    return dotprod;
}

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::norm(const Vector& vec) {
    T normret = 0.0;
    for (const auto& x:vec) {
        normret += x * x;
    }
    return sqrt(normret);
}

template<typename Vector, typename T>
inline T Lanczos<Vector, T>::l2norm(const Vector& alpha, const Vector& beta) {
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

template<typename Vector, typename T>
inline Vector& Lanczos<Vector, T>::normalise(Vector& vec) {
    T normal = norm(vec);
    for (auto& x:vec) {
        x /= normal;
    }
    return vec;
}

template<typename Vector, typename T>
Vector Lanczos<Vector, T>::init(const int& size) {
    Vector vec(size);
    for (auto& x:vec) {
        x = drand48();
    }
    T normalise = norm(vec);
    for (auto& x:vec) {
        x /= normalise;
    }

    //vec = {0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5};
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
