/*
 * =====================================================================================
 *
 *       Filename:  graph.cc
 *
 *    Description:  Member functions for Class Graph.hpp
 *        Created:  06/29/2016 18:42:42
 * 
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <exception>
#include <random>
#include <chrono>
#include <climits>

#include "graph.h"

using namespace std;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Constructor
 *  Description:  Generate random graph using addEdge function
 * =====================================================================================
 */

Graph::Graph(int num_of_vertex) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine rng(seed);
	uniform_int_distribution<int> num_of_neigh(1, 3);

	for (int vertex = 0; vertex < num_of_vertex; vertex++) {
		int num_of_neighbour = num_of_neigh(rng);
		//cout << "num_of_neighbour = " << num_of_neighbour << endl;
		uniform_int_distribution<int> randneigh(0, num_of_vertex - 1);
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
 *         Name:  init
 *  Description:  Initilise the subgraph
 * =====================================================================================
 */

void Graph::init(int rank, int global_size, int local_size) {
	rank_ = rank;
	global_size_ = global_size;
	local_size_ = local_size;
}

const int Graph::globalSize() const {
	return global_size_;
}

const int Graph::localSize() const {
	return local_size_;
}

const int Graph::rank() const {
	return rank_;
}

const int Graph::globalIndex(int local_index) const {
	return local_index + rank_ * local_size_;
}

const int Graph::localIndex(int global_index) const {
	return global_index - rank_ * local_size_;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  addEdge
 *  Description:  Add edges into graph 
 * =====================================================================================
 */

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
	/*	
	// Since graph is undirected, add edge dest->src	
	it = G.find(dest);
	if (it != G.end())
	it->second.insert(src);
	else {
	unordered_set<int> edges;
	edges.insert(src);
	G.insert(make_pair(dest, edges));
	}
	*/
}	

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  return graph properties
 *  Description:  Functions to return graph properties
 * =====================================================================================
 */

const int Graph::size() const {
	return G.size();
}

const int Graph::edgesNum() const {
	int edges = 0;
	for (auto& it:G) {
		edges += it.second.size(); 
	}
	return edges/2;
}

const int Graph::subgraphsNum() const {
	if (Colour.size() == 0) {
		cout << "WARNING:The graph hasn't been partitioned." << endl; 
		return 1;
	}

	unordered_map<int, int> reverse_colour_map;
	for (const auto& it:Colour) {
		reverse_colour_map.insert({it.second, it.first});
	}
	return reverse_colour_map.size();
}

const unordered_map<int, std::unordered_set<int>>::const_iterator Graph::find(int vertex) const {
	return G.find(vertex);
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
			cout << globalIndex(vertex) << ";" << endl;
		}
	else
		for (int vertex = 0; vertex < num_of_vertex; vertex++) {
			cout << globalIndex(vertex) << "[Colour=" << getColour(globalIndex(vertex)) << "];" << endl;
		}

	for (int vertex = 0; vertex < num_of_vertex; vertex++) {
		auto it = G.find(globalIndex(vertex));
		for (const int& neighbour:it->second)
			cout << globalIndex(vertex) << "--" << neighbour << " ;" << endl;
	}
	cout << "}" << endl;
}

/*-----------------------------------------------------------------------------
 *  Modified function to output the graph into dot file
 *-----------------------------------------------------------------------------*/

void Graph::outputDotFormat(const string& filename) const {
	int num_of_vertex = G.size();
	ofstream Output(filename);

	Output << "Undirected Graph {" << endl;
	if (Colour.size() == 0)
		for (int vertex = 0; vertex < num_of_vertex; vertex++) {
			Output << globalIndex(vertex) << ";" << endl;
		}
	else
		for (int vertex = 0; vertex < num_of_vertex; vertex++) {
			Output << globalIndex(vertex) << "[Colour=" << getColour(globalIndex(vertex)) << "];" << endl;
		}

	for (int vertex = 0; vertex < num_of_vertex; vertex++) {
		auto it = G.find(globalIndex(vertex));
		for (const int& neighbour:it->second)
			Output << globalIndex(vertex) << "--" << neighbour << " ;" << endl;
	}
	Output << "}" << endl;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printLaplacianMat
 *  Description:  Print the graph in Laplacian Matrix
 * =====================================================================================
 */

void Graph::printLaplacianMat() const {
	int num_of_vertex = G.size();

	int start = rank_ * num_of_vertex;
	int end = start + num_of_vertex;

	cout << "Laplacian Matrix:" << endl;
	for (int vertex = 0; vertex < end; vertex++) {
		cout << "\t";
		if (vertex >= start)
			cout << vertex;
	}
	cout << endl;

	for (int row = 0; row < num_of_vertex; row++) {
		cout << globalIndex(row) << "\t";
		auto it = G.find(globalIndex(row));
		for (int col = 0; col < global_size_; col++) {
			if (col == globalIndex(row))
				cout << G.at(globalIndex(row)).size() << "\t";
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
 *         Name:  setColour
 *  Description:  Map colour for each vertex
 * =====================================================================================
 */

void Graph::setColour(int vertex, int colour) const {
	//Colour.insert({vertex, colour});
	Colour[vertex] = colour;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getColour
 *  Description:  Return the colour of the vertex
 * =====================================================================================
 */

const int Graph::getColour(int vertex) const {
	if (Colour.size() == 0) {
		return 0;
	}
	return Colour.at(vertex);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  readDotFormat
 *  Description:  Read the graph from Dot file
 * =====================================================================================
 */

void Graph::readDotFormat(const string& filename) {

	ifstream In(filename);
	if (!In.is_open()) {
		std::cerr << "ERROR: Can't open the file" << endl;
		exit(-1);
	}
	In.ignore(INT_MAX, '-');
	In.ignore(1); // Skip the second '-'
	int from = 0, to = 0;
	In >> to;
	In.ignore(INT_MAX, '\n'); // Ignore other chars before end of line, go to next line

	while (In.good()) {
		if (from >= rank_ * local_size_ && from <= (rank_ + 1) * local_size_ - 1) {
			//cout << "from = " << from << " to = " << to << endl;
			addEdge(from, to);
			In >> from;
			In.ignore(2); // Ignore "--"
			In >> to;
			In.ignore(INT_MAX, '\n');
		}
		else {
			In >> from;
			In.ignore(2); // Ignore "--"
			In >> to;
			In.ignore(INT_MAX, '\n');
			if (from > (rank_ + 1) * local_size_ - 1) break;
		}
	}
	In.close();
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  readDotFormatWithColour
 *  Description:  Read the graph from Dot file where has colour for each vertex
 * =====================================================================================
 */

void Graph::readDotFormatWithColour(const string& filename) {

	ifstream In(filename);
	if (!In.is_open()) {
		std::cerr << "ERROR: Can't open the file" << endl;
		exit(-1);
	}
	In.ignore(INT_MAX, '{'); // Ignore the chars before the value of colour
	int vertex = 0, colour = 0;
	In >> vertex;
	int first_vertex = vertex;
	cout << "vertex = " << vertex << endl;
	In.ignore(INT_MAX, '='); // Ignore the chars before the value of colour
	In >> colour;
	cout << "colour = " << colour << endl;
	In.ignore(INT_MAX, '\n'); // Ignore other chars before end of line, go to next line

	local_size_ = G.size();
	while (In.good()) {
		setColour(vertex, colour);
		local_size_++;
		In >> vertex;
		if (vertex == first_vertex) break; // Break the loop at the beginning of edge line
		In.ignore(INT_MAX, '='); 
		In >> colour;
		In.ignore(INT_MAX, '\n');
	}

	int from = first_vertex, to = 0;
	cout << "first vertex = " << first_vertex << endl;
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
	In.close();
}
