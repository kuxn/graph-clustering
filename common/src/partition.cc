/**
 * @file partition.cc
 * @brief Partition the graph according to the eigenvectors
 * @author Ken Hu, xnchnhu@gmail.com
 */

#include "partition.h"
#include "lanczos.h"
#include "tqli.h"

#include <algorithm>
#include <boost/timer.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#ifdef VT_
#include "vt_user.h"
#endif

#define Sign(a) (a >= 0.0 ? 1 : 0)

using namespace std;

using DenseMatrix = std::vector<std::vector<double>>;

/**
 * @brief Partition the graph into multiple subgraphs by only calculating the
 *        corresponding laplacian eigenvectors
 * @param g The graph to partition
 * @param numOfSubGraphs The number of subgraphs to partition to
 * @param enableGramSchmidt Enable GramSchmidt
 */
Partition::Partition(const Graph& g, const int& numOfSubGraphs, bool enableGramSchmidt)
{
#ifdef VT_
    VT_TRACER("Partition::Partition");
#endif
    boost::timer timer_partition;
    int numOfEigenvectors = log2(numOfSubGraphs);

    // Construct tridiagonal matrix using Lanczos algorithm
    boost::timer lanczosTimer;
    Lanczos<std::vector<double>, double> lanczos(g, numOfEigenvectors, enableGramSchmidt);
    double t_lan = lanczosTimer.elapsed();
    times.push_back(t_lan);
    laplacianEigenvalues_ = lanczos.alpha;
    std::vector<double> beta = lanczos.beta;

    // Define an identity matrix as the input for TQLI algorithm
    DenseMatrix tridiagonalEigenvectors;

    // Calculate the eigenvalues and eigenvectors of the tridiagonal matrix
    boost::timer tqliTimer;
    tqli(laplacianEigenvalues_, beta, tridiagonalEigenvectors);
    double t_tqli = tqliTimer.elapsed();
    times.push_back(t_tqli);

    // Find the index of the nth smallest eigenvalue (fiedler vector) of the
    // eigenvalues vector "alpha"
    int vectorIndex = 0;

    int m = laplacianEigenvalues_.size();
    unordered_multimap<double, int> hashmap;
    for (int i = 0; i < m; i++) {
        hashmap.insert({laplacianEigenvalues_[i], i});
    }
    std::vector<double> auxiliaryVector = laplacianEigenvalues_;
    sort(auxiliaryVector.begin(), auxiliaryVector.end());

#ifndef Median_
    int fielderIndex = 1;
    for (int i = 0; i < numOfEigenvectors; i++) {
        auto it = hashmap.find(auxiliaryVector[fielderIndex]);
        while (abs(it->first) < 1e-2) {
            fielderIndex++;
            it = hashmap.find(auxiliaryVector[fielderIndex]);
        }
        fielderIndex++;
        vectorIndex = it->second;
        ritzValues.push_back(it->first);
        hashmap.erase(it);  // Deal with identical eigenvalues
        laplacianEigenMatrix_.push_back(getOneLapEigenVec(
            lanczos.lanczos_vecs, tridiagonalEigenvectors, vectorIndex));
    }

    for (int vertex = 0; vertex < g.size(); vertex++) {
        int colour = 0;
        for (int row = 0; row < numOfEigenvectors; row++) {
            colour += pow(2, row) * Sign(laplacianEigenMatrix_[row][vertex]);
        }
        g.setColour(g.globalIndex(vertex), colour);
    }
#endif

#ifdef Median_
    std::vector<double> medianVector;
    double median = 0.0;
    int fielderIndex = 1;
    for (int i = 0; i < numOfEigenvectors; i++) {
        auto it = hashmap.find(auxiliaryVector[fielderIndex]);
        while (abs(it->first) < 1e-2) {
            fielderIndex++;
            it = hashmap.find(auxiliaryVector[fielderIndex]);
        }
        fielderIndex++;
        vectorIndex = it->second;
        ritzValues.push_back(it->first);
        hashmap.erase(it);  // Deal with identical eigenvalues
        laplacianEigenMatrix_.push_back(getOneLapEigenVec(
            lanczos.lanczos_vecs, tridiagonalEigenvectors, vectorIndex));

        // Calculate the median for each eigenvector
        std::vector<double> auxiliaryVector2 = laplacianEigenMatrix_[i];
        sort(auxiliaryVector2.begin(), auxiliaryVector2.end());
        int vectorSize = auxiliaryVector2.size();
        if (vectorSize % 2 == 0) {
            median = (auxiliaryVector2[vectorSize / 2] +
                      auxiliaryVector2[vectorSize / 2 - 1]) /
                     2;
        } else {
            median = auxiliaryVector2[vectorSize / 2];
        }
        medianVector.push_back(median);
    }
    for (int vertex = 0; vertex < g.size(); vertex++) {
        int colour = 0;
        for (int row = 0; row < numOfEigenvectors; row++) {
            colour +=
                pow(2, row) *
                signMedian(laplacianEigenMatrix_[row][vertex], medianVector[row]);
        }
        g.setColour(g.globalIndex(vertex), colour);
    }
#endif

    double t_par = timer_partition.elapsed();
    times.push_back(t_par);
}

inline int Partition::signMedian(double entry, double median)
{
    return entry >= median ? 1 : 0;
}

/**
 * @brief Print eigenvalues, eigenvectors
 */
void Partition::printLapEigenMat()
{
    int row_size = laplacianEigenMatrix_.size();
    for (int row = 0; row < row_size; row++) {
        int col_size = laplacianEigenMatrix_[row].size();
        for (int col = 0; col < col_size; col++) {
            cout << laplacianEigenMatrix_[row][col] << " ";
        }
        cout << endl;
    }
}

void Partition::printLapEigenvalues()
{
    cout << "Laplacian Eigenvalues: ";
    for (const double& x : laplacianEigenvalues_) {
        cout << x << " ";
    }
    cout << endl;
}

void Partition::outputLapEigenvalues()
{
    string filename("./output/eigenvalues_");
    filename += to_string(laplacianEigenvalues_.size());
    filename += ".dat";

    ofstream Output(filename);

    for (unsigned int i = 0; i < laplacianEigenvalues_.size(); i++) {
        Output << i << " " << laplacianEigenvalues_[i] << endl;
    }

    Output.close();
}

/**
 * @brief Calculate one eigenvector.
 * @param lanczosVectors
 * @param tridiagonalEigenvectors
 * @param vectorIndex
 * @return laplacianVectors
 */
std::vector<double> Partition::getOneLapEigenVec(DenseMatrix& lanczosVectors,
                                    DenseMatrix& tridiagonalEigenvectors,
                                    const int& vectorIndex)
{
#ifdef VT_
    VT_TRACER("Partition::getOneLapEigenVec");
#endif
    // Find the eigenvector from the vector matrix of the Tridiagonal
    // eigenvector matrix (each column represents a vector)
    int size = tridiagonalEigenvectors.size();  // Size of vertices/vectors
    std::vector<double> vec(size, 0);
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            if (col == vectorIndex) {
                vec[row] = tridiagonalEigenvectors[row][col];
            }
        }
    }
    // Calculate the corresponding Laplacian vector by Lanczos vectors(each row
    // represents a vector, need transposing)
    // lanczos vector - m * n, tridiagonalEigenvectors - m * m, Ritz - m * n, Ritz_vector
    // - n
    int col_size = lanczosVectors[0].size();
    int row_size = lanczosVectors.size();
    std::vector<double> laplacianVector(col_size, 0);
    for (int col = 0; col < col_size; col++) {
        for (int row = 0; row < row_size; row++) {
            laplacianVector[col] +=
                lanczosVectors[row][col] * vec[row];
        }
    }
    return laplacianVector;
}
