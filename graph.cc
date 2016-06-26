/*
 * =====================================================================================
 *
 *       Filename:  graph.cpp
 *
 *    Description:  Member functions for Class Graph.hpp
 *        Created:  06/09/2016 22:42:42
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

const int Graph::edgesNum() const {
	return edges_/2;
}

const unsigned int Graph::size() const {
	return G.size();
}

unordered_map<int, std::unordered_set<int>>::const_iterator Graph::find(int vertex) const {
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
    int num_of_vertex = G.size();
	cout << "Undirected Graph {" << endl;
	if (Colour.size() == 0)
		for (int vertex = 0; vertex < num_of_vertex; vertex++) {
				cout << vertex << ";" << endl;
		}
	else
		for (int vertex = 0; vertex < num_of_vertex; vertex++) {
			cout << vertex << "[Colour=" << getColour(vertex) << "];" << endl;
		}

	for (int vertex = 0; vertex < num_of_vertex; vertex++) {
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
	int num_of_vertex = G.size();
	cout << "Laplacian Matrix:" << endl;
	for (int vertex = 0; vertex < num_of_vertex; vertex++) {
		cout << "\t" << vertex;
	}
	cout << endl;
	
	for (int row = 0; row < num_of_vertex; row++) {
		cout << row << "\t";
		auto it = G.find(row);
		for (int col = 0; col < num_of_vertex; col++) {
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

Graph::Graph(int num_of_vertex) {
	random_device seed;	mt19937 rng (seed());
	//uniform_int_distribution<int> num_of_neigh(0, num_of_vertex - 1);
	uniform_int_distribution<int> num_of_neigh(1, 3);
	//poisson_distribution<int> num_of_neigh(num_of_vertex/2);

	edges_ = 0;
	for (int vertex = 0; vertex < num_of_vertex; vertex++) {
		int num_of_neighbour = num_of_neigh(rng);
		edges_ += num_of_neighbour;
        //cout << "num_of_neighbour = " << num_of_neighbour << endl;
		uniform_int_distribution<int> randneigh(0, num_of_vertex - 1);
		//poisson_distribution<int> randneigh(num_of_vertex/2);
		unordered_set<int> neighbours;
		for (int neighbour = 0; neighbour < num_of_neighbour; neighbour++) {
			int rand_neighbour = 0;
			int trials = 1000;
			do {
				rand_neighbour = randneigh(rng);
				trials--;
			} while (vertex == rand_neighbour && trials);
			
			if (trials <= 0) throw std::runtime_error ("Run out of trials.");
			addEdge(vertex, rand_neighbour);
		}
	}
	if (G.size() != (unsigned)num_of_vertex) throw std::length_error ("The size of generated graph is incorrect.");
	cout << "Graph generation is done." << endl;
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
	edges_ = 0;
	In >> to;
	In.ignore(INT_MAX, '\n'); // Ignore other chars before end of line, go to next line
	while (In.good()) {
		addEdge(from, to);
		edges_++;
		In >> from;
		In.ignore(2); // Ignore "--"
		In >> to;
		In.ignore(INT_MAX, '\n');
	}
	In.close();
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
	edges_ = 0;
	In.ignore(2); // Ignore "--"
	In >> to;
	In.ignore(INT_MAX, '\n'); 
	while (In.good()) {
		addEdge(from, to);
		edges_++;
		In >> from;
		In.ignore(2); // Ignore "--"
		In >> to;
		In.ignore(INT_MAX, '\n');
	}
	In.close();
}
