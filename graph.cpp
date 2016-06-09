/*
 * =====================================================================================
 *
 *       Filename:  graph.cpp
 *
 *    Description:  Member functions for Class Graph.hpp
 *
 * 		  Version:  1.0
 *		  Created:  06/09/2016 22:42:42
 *		 Revision:  none
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
	else
	{
		setOfEdges edges;
		edges.insert(dest);
		G.insert(make_pair(src, edges));
	}
	
	// Since graph is undirected, add edge dest->src	
	it = G.find(dest);
	if (it != G.end())
		it->second.insert(src);
	else
	{
		setOfEdges edges;
		edges.insert(src);
		G.insert(make_pair(dest, edges));
	}
	
}	


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  print
 *  Description:  Print graph in DOT format
 * =====================================================================================
 */

void Graph::print() {
	cout << "Undirected Graph {" << endl;
	for (unsigned int indexOfVertex = 0; indexOfVertex < G.size(); indexOfVertex++) {
		cout << indexOfVertex << ";" << endl;
	}

	for (unsigned int indexOfVertex = 0; indexOfVertex < G.size(); indexOfVertex++) {
		it = G.find(indexOfVertex);
		for (const int& vertex:it->second)
			cout << indexOfVertex << "--" << vertex << " ;" << endl;
	}
	cout << "}" << endl;
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
	uniform_int_distribution<int> randEdgeSet(0, numOfVertex - 1);

	for (int indexOfVertex = 0; indexOfVertex < numOfVertex; indexOfVertex++) {
		int sizeOfEdgeSet = randEdgeSet(rng);
        cout << "sizeOfEdgeSet = " << sizeOfEdgeSet << endl;
		uniform_int_distribution<int> randEdge(0, sizeOfEdgeSet);
		setOfEdges edges;
		for (int indexOfEdge = 0; indexOfEdge < sizeOfEdgeSet; indexOfEdge++) {
			int randEdgeIndex = randEdge(rng);
			if (randEdgeIndex != indexOfVertex)
			    addEdge(indexOfVertex, randEdgeIndex);
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
		setOfEdges edges;
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
