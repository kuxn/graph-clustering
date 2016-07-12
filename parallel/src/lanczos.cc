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
 *  Modified Lanczos algorithm with GramSchmidt by Gram–Schmidt
 *-----------------------------------------------------------------------------*/

template<typename Vector, typename T>
Lanczos<Vector, T>::Lanczos(const Graph& g_local, bool GramSchmidt) {

	VT_TRACER("LANCZOS");
	int local_size = g_local.localSize();
	Vector v0_local(local_size);
	Vector v1_halo(g_local.globalSize());

	v0_local = init(v0_local, g_local);

#ifdef Debug
	cout << "rank " << g_local.rank() << ": ";
	for (int i = 0; i < local_size; i++) {
		cout << "v0_local[" << i << "] = " << v0_local[i] << " ";
	}
	cout << endl;
	cout << "global dot of these vectors = " << dot(v0_local, v0_local) << endl;
#endif

	Vector t_local = v0_local, v1_local = v0_local, w_local;
	T alpha_val_global = 0.0, beta_val_global = 0.0;

#ifdef Debug
	cout << "IN rank " << g_local.rank() << " before halo update v1_halo:";
	for (int i = 0; i < g_local.globalSize(); i++) {
		cout << v1_halo[i] << " ";
	}
	cout << endl;
#endif

	lanczos_vecs_local[0] = v1_local;
	haloInit(g_local);

#ifdef Debug
	cout << "IN rank " << g_local.rank() << " before halo update v1_halo:";
	cout << "IN rank " << g_local.rank() << " after halo update v1_halo:";
	for (int i = 0; i < g_local.globalSize(); i++) {
		cout << v1_halo[i] << " ";
	}
	cout << endl;
#endif

	for (int iter = 1; iter < g_local.globalSize(); iter++) {
		haloUpdate(g_local, v1_local, v1_halo);

#ifdef Debug
		if (g_local.rank() == 0) {
			cout << "IN rank " << g_local.rank() << " after " << iter <<" halo update v1_halo:";
			for (int i = 0; i < g_local.globalSize(); i++) {
				cout << v1_halo[i] << " ";
			}
			cout << endl;
		}
#endif
		w_local = multGraphVec(g_local, v1_halo);

#ifdef Debug
		cout << "w_local, rank " << g_local.rank() << ": ";
		for (auto x:w_local)
			cout << x << " ";
		cout << endl;
#endif

		T alpha_val_global = dot(v1_local, w_local);
		//cout << "dot(v1, v1) = " << dot(v1, v1) << endl;
		alpha_global.push_back(alpha_val_global);

#ifdef Debug
		cout << "rank " << g_local.rank() << ": ";
		cout << "alpha[" << iter - 1 << "] =  " << alpha_val_global << endl;
#endif

		for (int i = 0; i < local_size; i++) {
			t_local[i] = w_local[i] - alpha_val_global * v1_local[i] - beta_val_global * v0_local[i];
		}

		beta_val_global = sqrt(dot(t_local, t_local)); 
		beta_global.push_back(beta_val_global);	
		//cout << "beta[" << iter-1 << "]: " << beta_val_global << endl;

		if (std::abs(beta_val_global) < 1e-5) 
			try { throw std::runtime_error("Value of beta is close to 0: "); }
		catch (std::runtime_error& e) { 
			std::cerr << "ERROR: " << e.what(); 
			cout << "beta[" << iter-1 << "]: " << beta_val_global << endl;
		}

		if (!GramSchmidt) {
			v0_local = v1_local;
		}
		for (int i = 0; i < local_size; i++) {
			v1_local[i] = t_local[i]/beta_val_global;
		}

		lanczos_vecs_local[iter] = v1_local;

		if (GramSchmidt) {
			boost::mpi::timer timer_gram_schmidt;
			gramSchmidt(iter, g_local);
			v0_local = lanczos_vecs_local[iter-1];
			v1_local = lanczos_vecs_local[iter];
		}

		//Verify the dot product of v0 and v1 which is supposed to be 0
		T dot_global = dot(v0_local, v1_local);
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
	haloUpdate(g_local, v1_local, v1_halo);
	w_local = multGraphVec(g_local, v1_halo);
	alpha_val_global = dot(v1_local, w_local);
	alpha_global.push_back(alpha_val_global);

	transform(g_local);	
	if (GramSchmidt && g_local.rank() == 0) {
		cout << "Lanczos algorithm WITH GramSchmidt is done." << endl;
	} else if (g_local.rank() == 0){
		cout << "Lanczos algorithm WITHOUT GramSchmidt is done." << endl; 
	}
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
#ifdef Debug
	for (int rank = 0; rank < world.size(); rank++) {
		if (rank != g.rank()) {
			auto it = halo_recv.find(rank);
			cout <<"rank " << g.rank() << " receive from rank " << rank << " - "; 
			for (const int& halo_neighbour:it->second) {
				cout << halo_neighbour << " ";
			}
			cout << endl;
		}
	}
	for (int rank = 0; rank < world.size(); rank++) {
		if (rank != g.rank()) {
			auto it = halo_send.find(rank);
			cout <<"rank " << g.rank() << " send to rank " << rank << " - "; 
			for (const int& halo_neighbour:it->second) {
				cout << halo_neighbour << " ";
			}
			cout << endl;
		}
	}
#endif
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  haloUpdate
 *  Description:  Refresh the halo elements each iteration for Graph * Lanczos_Vec
 * =====================================================================================
 */

template<typename Vector, typename T>
void Lanczos<Vector, T>::haloUpdate(const Graph& g, Vector& v_local, Vector& v_halo) {
	mpi::request reqs[world.size()-1];
	int i = 0;
	for (int rank = 0; rank < world.size(); rank++) {
		if (rank != g.rank()) {
			auto it = halo_send.find(rank);
			if (it != halo_send.end()) {
				std::unordered_map<int, T> buf_send(it->second.size()); // <global_index, value>;
				for (const int& halo_neighbour:it->second) {
					buf_send.insert({halo_neighbour, v_local[g.localIndex(halo_neighbour)]}); // global index, locale value
				}
				reqs[i] = world.isend(rank, 0, buf_send); //(dest, tag, value to send)
				i++;
			}
		} 
	}
	std::vector<std::unordered_map<int, T>> buf_recv(world.size()); // <global_index, value>;
	//std::unordered_map<int, std::unordered_map<int, T>> buf_recv; // <global_index, value>;
	i = 0;
	for (int rank = 0; rank < world.size(); rank++) {
		if (rank != g.rank()) {
			auto it = halo_recv.find(rank);
			if (it != halo_recv.end()) {
				reqs[i] = world.irecv(rank, 0, buf_recv[i]); //(src, tag, store to value)
				//cout << "in rank " << g.rank() << " buf["<< rank << "] received from rank " << rank << endl; 
				i++;
			}
		} else {
			for (int j = 0; j < g.localSize(); j++) {
				v_halo[g.globalIndex(j)] = v_local[j];	
			}
		}
	}
	mpi::wait_all(reqs, reqs+world.size()-1);
	std::unordered_map<int, T> buf_recv_one_pack; // <global_index, value>;
	i = 0;
	// Unpack the buffer to fill in to v_halo
	for (int rank = 0; rank < world.size(); rank++) {
		if (rank != g.rank()) {
			//cout << "in rank " << g.rank() << " buf["<< rank << "] unpaking data from rank " << rank << endl; 
			buf_recv_one_pack = buf_recv[i];
			auto it = halo_recv.find(rank);
			if (it != halo_recv.end()) {
				for (const int& halo_neighbour:it->second) {
					v_halo[halo_neighbour] = buf_recv_one_pack[halo_neighbour]; //(src, tag, store to value)
				}
			i++;
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
inline void Lanczos<Vector, T>::gramSchmidt(int& iter, const Graph& g) {
    //VT_TRACER("GramSchmidt");
    for (int k = 1; k <= iter; k++) {
        //cout << "i - norm of lanczos_vecs_local["<<k<<"] = " << norm(lanczos_vecs_local[k]) << endl;
        for (int i = 0; i < k; i++) {
            T dot_global = dot(lanczos_vecs_local[i], lanczos_vecs_local[k]);
            //if (world.rank() == 0) {
            //	cout << "iter " << iter << " gramSchmidt global dot product " << dot_global << endl;
            //}
            for (int j = 0; j < g.localSize(); j++) {
                lanczos_vecs_local[k][j] -= dot_global * lanczos_vecs_local[i][j];
            }
        }
        T norm_global = std::sqrt(dot(lanczos_vecs_local[k], lanczos_vecs_local[k]));
        //if (g.rank() == 0 && (norm_global - 1) > 1e-2) {
        //    cout << "global dot of lanczos_vecs_local[k] = " << norm_global<< endl;
        //}
        normalise(lanczos_vecs_local[k], g, norm_global);
        //cout << "local dot of lanczos_vecs_local[k] = " << dot_local(lanczos_vecs_local[k], lanczos_vecs_local[k]) << endl;
        //cout << "global dot of lanczos_vecs_local[k] = " << dot(world, lanczos_vecs_local[k], lanczos_vecs_local[k]) << endl;
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
	for (const T& value:vec) {
		norm_local += value * value;
	}
	return sqrt(norm_local);
}

template<typename Vector, typename T>
inline void Lanczos<Vector, T>::normalise(Vector& vec, const Graph& g, const T& norm_global) {
	int local_size = vec.size();
	for (int i = 0; i < local_size; i++) {
		vec[i] /= norm_global;
	}
}

template<typename Vector, typename T>
Vector& Lanczos<Vector, T>::init(Vector& vec, const Graph& g) {
	int local_size = vec.size();
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> gen(0.0,1.0);
	for (int i = 0; i < local_size; i++) {
		vec[i] = gen(generator);
		//vec[i] = i;
	}
	T norm_local = norm(vec);
	for (int i = 0; i < local_size; i++) {
		vec[i] /= norm_local;
		vec[i] /= sqrt(g.globalSize()/local_size); // vec[i]/=sqrt(procs), to make sure the global vector is normalised
	}
	return vec;
}

template<typename Vector, typename T>
void Lanczos<Vector, T>::transform(const Graph& g) {
#ifdef Debug
	for (unsigned int i = 0; i < lanczos_vecs_local.size(); i++) {
		auto it = lanczos_vecs_local.find(i);
		cout << "rank " << g.rank() << " ";
		for (auto x:it->second) {
			cout << x << "\t";
		}
		cout << endl;
	}
#endif
	std::vector<std::unordered_map<int, Vector>> gather;
	all_gather(world, lanczos_vecs_local, gather);
	//cout << "size of all_gather = " << gather.size() << endl;

	Vector vinitial(g.globalSize(), 0);
	for(int i = 0; i < g.globalSize(); i++)	lanczos_vecs_global[i] = vinitial;

	for (int row = 0; row < g.globalSize(); row ++) {
		for (int i = 0; i < world.size(); i++) {
			for (int j = 0; j < g.localSize(); j++) {
#ifdef Debug
				if (g.rank() == 0) {
					cout << "lanczos_vecs_global[" << row << "][" << i  * g.localSize() + j << "] = "; 
					cout << "gather[" << i << "][" << j << "][" << row << "] = " << gather[i][row][j] << endl; 
				}
#endif
				lanczos_vecs_global[row][i * g.localSize() + j] = gather[i][row][j];
			}
		}
	}
#ifdef Debug
	if (g.rank() == 0) {
		for (unsigned int i = 0; i < lanczos_vecs_global.size(); i++) {
			auto it = lanczos_vecs_global.find(i);
			cout << "rank " << g.rank() << " ";
			for (auto x:it->second) {
				cout << x << "\t";
			}
			cout << endl;
		}
	}
#endif
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
