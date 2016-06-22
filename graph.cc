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
#include <stdexcept>
#include <random>
#include <climits>

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

	// Avoid the self circle
	if (src == dest) return;

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
	if (Colour.size() == 0)
		for (int vertex = 0; vertex < numofvertex; vertex++) {
				cout << vertex << ";" << endl;
		}
	else
		for (int vertex = 0; vertex < numofvertex; vertex++) {
			cout << vertex << "[Colour=" << getColour(vertex) << "];" << endl;
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
	//uniform_int_distribution<int> numofneigh(0, numofvertex - 1);
	uniform_int_distribution<int> numofneigh(1, 3);
	//poisson_distribution<int> numofneigh(numofvertex/2);

	for (int vertex = 0; vertex < numofvertex; vertex++) {
		int numofneighbour = numofneigh(rng);
        //cout << "numofneighbour = " << numofneighbour << endl;
		uniform_int_distribution<int> randneigh(0, numofvertex - 1);
		//poisson_distribution<int> randneigh(numofvertex/2);
		unordered_set<int> neighbours;
		for (int neighbour = 0; neighbour < numofneighbour; neighbour++) {
			int randneighbour = 0;
			do {
				randneighbour = randneigh(rng);
			} while (vertex == randneighbour);
			addEdge(vertex, randneighbour);
		}
	}
	if (G.size() != (unsigned)numofvertex) throw std::length_error ("The size of generated graph is incorrect.");
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  setColour
 *  Description:  Map colour for each vertex
 * =====================================================================================
 */

void Graph::setColour(int vertex, int colour) const {
	Colour.insert({vertex, colour});
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getColour
 *  Description:  Return the colour of the vertex
 * =====================================================================================
 */

const int Graph::getColour(int vertex) const {
	return Colour.at(vertex);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  readDotFormat
 *  Description:  Read the graph from Dot file
 * =====================================================================================
 */

void Graph::readDotFormat(ifstream& In) {
	In.ignore(INT_MAX, '-');
	In.ignore(1); // Skip the second '-'

	int from = 0;
	int to = 0;
	In >> to;
	In.ignore(INT_MAX, '\n'); // Ignore other chars before end of line, go to next line
	while (In.good()) {
		addEdge(from, to);
		In >> from;
		In.ignore(2); // Ignore "--"
		In >> to;
		In.ignore(INT_MAX, '\n');
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  readDotFormatWithColour
 *  Description:  Read the graph from Dot file where has colour for each vertex
 * =====================================================================================
 */

void Graph::readDotFormatWithColour(ifstream& In) {
	
	In.ignore(INT_MAX, '='); // Ignore the chars before the value of colour
	int vertex = 0;
	int colour = 0;
	In >> colour;
	In.ignore(INT_MAX, '\n'); // Ignore other chars before end of line, go to next line
	
	while (In.good()) {
		setColour(vertex, colour);
		In >> vertex;
		if (vertex == 0) break; // Break the loop at the beginning of edge line
		In.ignore(INT_MAX, '='); 
		In >> colour;
		In.ignore(INT_MAX, '\n');
	}

	int from = 0;
	int to = 0;
	In.ignore(2); // Ignore "--"
	In >> to;
	In.ignore(INT_MAX, '\n'); 
	while (In.good()) {
		addEdge(from, to);
		In >> from;
		In.ignore(2); // Ignore "--"
		In >> to;
		In.ignore(INT_MAX, '\n');
	}
}
