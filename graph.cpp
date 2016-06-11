/*
 * =====================================================================================
 *
 *       Filename:  graph.cpp
 *
 *    Description:  Member functions for Class Graph.hpp
 *
 * 		  Version:  1.1
 *		  Created:  06/09/2016 22:42:42
 *		 Revision:  Add printLaplacianMat()
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <random>

#include "graph.hpp"

using namespace std;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  addEdge
 *  Description:  Add edges into graph 
 * =====================================================================================
 */

void Graph::addEdge(int src, int dest) {

	// Add edge src->edge
	it = G.find(src);
	if (it != G.end())
		it->second.insert(dest);
	else {
		unordered_set<int> edges;
		edges.insert(dest);
		G.insert(make_pair(src, edges));
	}
	
	// Since graph is undirected, add edge dest->src	
	it = G.find(dest);
	if (it != G.end())
		it->second.insert(src);
	else {
		unordered_set<int> edges;
		edges.insert(src);
		G.insert(make_pair(dest, edges));
	}
	
}	

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printDotFormat
 *  Description:  Print graph in DOT format
 * =====================================================================================
 */

void Graph::printDotFormat() {
    numOfVertex = G.size();
	cout << "Undirected Graph {" << endl;
	for (int indexOfVertex = 0; indexOfVertex < numOfVertex; indexOfVertex++) {
		cout << indexOfVertex << ";" << endl;
	}

	for (int indexOfVertex = 0; indexOfVertex < numOfVertex; indexOfVertex++) {
		it = G.find(indexOfVertex);
		for (const int& neighbour:it->second)
			cout << indexOfVertex << "--" << neighbour << " ;" << endl;
	}
	cout << "}" << endl;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printLaplacianMat
 *  Description:  Print the graph in Laplacian Matrix
 * =====================================================================================
 */

void Graph::printLaplacianMat() {
    numOfVertex = G.size();
    cout << "Laplacian Matrix:" << endl;
	for (int indexOfVertex = 0; indexOfVertex < numOfVertex; indexOfVertex++) {
        cout << "\t" << indexOfVertex;
    }
    cout << endl;

	for (int indexOfRow = 0; indexOfRow < numOfVertex; indexOfRow++) {
        cout << indexOfRow << "\t";
        it = G.find(indexOfRow);
        for (int indexOfCol = 0; indexOfCol < numOfVertex; indexOfCol++) {
            if (indexOfCol == indexOfRow)
                cout << G.at(indexOfRow).size() << "\t";
            else if (it->second.find(indexOfCol) != it->second.end())
                cout << "-1\t";
            else
                cout << "0\t";
        }
        cout << endl;
    }    
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  genRandomGraph
 *  Description:  Generate random graph using addEdge function
 * =====================================================================================
 */

void Graph::genRandomGraph(int num) {
	random_device seed;	mt19937 rng (seed());
	numOfVertex = num;
	uniform_int_distribution<int> NumOfNeighbours(0, numOfVertex - 1);

	for (int indexOfVertex = 0; indexOfVertex < numOfVertex; indexOfVertex++) {
		int numOfNeighbours = NumOfNeighbours(rng);
        //cout << "sizeOfEdgeSet = " << sizeOfEdgeSet << endl;
		uniform_int_distribution<int> RandNeighbour(0, numOfNeighbours);
		unordered_set<int> neighbours;
		for (int indexOfNeighbour = 0; indexOfNeighbour < numOfNeighbours; indexOfNeighbour++) {
			int randNeighbour = RandNeighbour(rng);
			if (randNeighbour != indexOfVertex)
			    addEdge(indexOfVertex, randNeighbour);
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  genRandomGraph
 *  Description:  Generate random graph by generate random set for each vertex
 * =====================================================================================
 */

/*
void genRandomGraph(int num) {
	random_device seed;	mt19937 rng (seed());
	numOfVertex = num;
	//cout << "numOfVertex = " << numOfVertex << endl;
	for (int indexOfVertex = 0; indexOfVertex < numOfVertex; indexOfVertex++) {
		int sizeOfEdgeSet = rand() % numOfVertex;
		//cout << "sizeOfEdgeSet = " << sizeOfEdgeSet << endl;
		uniform_int_distribution<int> randEdge(0, sizeOfEdgeSet);
		unordered_map<int> edges;
		for (int indexOfEdge = 0; indexOfEdge < sizeOfEdgeSet; indexOfEdge++) {
			int randEdgeIndex = randEdge(rng);
			//cout << "randEdgeIndex = " << randEdgeIndex << endl;
			if (randEdgeIndex != indexOfVertex)
			edges.insert(randEdgeIndex);	
		}
		//cout << "indexOfVertex = " << indexOfVertex << endl;
		G[indexOfVertex] = edges;
	}
}
*/

