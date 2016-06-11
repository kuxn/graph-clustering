/*
 * =====================================================================================
 *
 *       Filename:  graph.hpp
 *
 *    Description:  Class Graph header file
 *
 * 		  Version:  1.1
 *		  Created:  06/09/2016 16:22:24
 *		 Revision:  Add printLaplacianMat
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

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <unordered_set>
#include <unordered_map>

using namespace std;

typedef unordered_set<int> SetOfNeighbours;

class Graph {

	private:
		int num_of_vertex = 10;
		unordered_map<int, SetOfNeighbours> G;

	public:
		Graph() {}
		Graph(int num): num_of_vertex(num) {}
		
		void addEdge(int src, int dest);
		void genRandomGraph(int num);
		void printDotFormat() const;
		void printLaplacianMat() const;
		const unsigned int size() const;
		unordered_map<int, SetOfNeighbours>::const_iterator find(int vertex) const;
};

#endif
