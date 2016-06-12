/*
 * =====================================================================================
 *
 *       Filename:  graph.cpp
 *
 *    Description:  Member functions for Class Graph.hpp
 * 	  	  Created:  06/09/2016 22:42:42
 * 
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <random>

#include "graph.h"

using namespace std;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  addEdge
 *  Description:  Add edges into graph 
 * =====================================================================================
 */

const unsigned int Graph::size() const {
	return G.size();
}

unordered_map<int, SetOfNeighbours>::const_iterator Graph::find(int vertex) const {
	return G.find(vertex);
}

void Graph::addEdge(int src, int dest) {

	// Add edge src->edge
	auto it = G.find(src);
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

void Graph::printDotFormat() const {
    int numofvertex = G.size();
	cout << "Undirected Graph {" << endl;
	for (int vertex = 0; vertex < numofvertex; vertex++) {
		cout << vertex << ";" << endl;
	}

	for (int vertex = 0; vertex < numofvertex; vertex++) {
		auto it = G.find(vertex);
		for (const int& neighbour:it->second)
			cout << vertex << "--" << neighbour << " ;" << endl;
	}
	cout << "}" << endl;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printLaplacianMat
 *  Description:  Print the graph in Laplacian Matrix
 * =====================================================================================
 */

void Graph::printLaplacianMat() const {
    int numofvertex = G.size();
    cout << "Laplacian Matrix:" << endl;
	for (int vertex = 0; vertex < numofvertex; vertex++) {
        cout << "\t" << vertex;
    }
    cout << endl;

	for (int row = 0; row < numofvertex; row++) {
        cout << row << "\t";
        auto it = G.find(row);
        for (int col = 0; col < numofvertex; col++) {
            if (col == row)
                cout << G.at(row).size() << "\t";
            else if (it->second.find(col) != it->second.end())
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

void Graph::genRandomGraph(int numofvertex) {
	random_device seed;	mt19937 rng (seed());
	uniform_int_distribution<int> numofneigh(0, numofvertex - 1);

	for (int vertex = 0; vertex < numofvertex; vertex++) {
		int numofneighbour = numofneigh(rng);
        //cout << "sizeofedgeset = " << sizeofedgeset << endl;
		uniform_int_distribution<int> randneigh(0, numofneighbour);
		unordered_set<int> neighbours;
		for (int neighbour = 0; neighbour < numofneighbour; neighbour++) {
			int randneighbour = randneigh(rng);
			if (randneighbour != vertex)
			    addEdge(vertex, randneighbour);
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
void genRandomGraph(int numofvertex) {
	random_device seed;	mt19937 rng (seed());
	//cout << "numofvertex = " << numofvertex << endl;
	for (int vertex = 0; vertex < numofvertex; vertex++) {
		int sizeofedgeset = rand() % numofvertex;
		//cout << "sizeofedgeset = " << sizeofedgeset << endl;
		uniform_int_distribution<int> randedge(0, sizeofedgeset);
		unordered_map<int> edges;
		for (int edge = 0; edge < sizeofedgeset; edge++) {
			int randedgeindex = randedge(rng);
			//cout << "randedgeindex = " << randedgeindex << endl;
			if (randedgeindex != vertex)
			edges.insert(randedgeindex);	
		}
		//cout << "vertex = " << vertex << endl;
		G[vertex] = edges;
	}
}
*/

