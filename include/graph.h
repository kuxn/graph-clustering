/*
 * =====================================================================================
 *
 *       Filename:  graph.h
 *
 *    Description:  Class Graph header file
 *        Created:  06/09/2016 16:22:24
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */


/*
 * =====================================================================================
 *        Class:  Graph
 *  Description:  Class to create a new graph object
 * =====================================================================================
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <unordered_set>
#include <unordered_map>
#include <fstream>

class Graph {

	private:
		int num_of_vertex = 10;
		typedef std::unordered_set<int> SetOfNeighbours;
		std::unordered_map<int, SetOfNeighbours> G;
		mutable std::unordered_map<int, int> Colour;

	public:
		Graph() {}
		Graph(int num): num_of_vertex(num) {}
		
		void addEdge(int src, int dest);
		void genRandomGraph(int num);
		void printDotFormat() const;
		void printLaplacianMat() const;
		const unsigned int size() const;
		std::unordered_map<int, std::unordered_set<int>>::const_iterator find(int vertex) const;
		void setColour(int vertex, int colour) const;
		const int getColour(int vertex) const;
		void readDotFormat(std::ifstream& In);
		void readDotFormatWithColour(std::ifstream& In);
};

#endif
