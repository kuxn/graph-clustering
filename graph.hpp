#ifndef GRAPH_HPP_
#define GRAPH_HPP_

/*
 * =====================================================================================
 *
 *       Filename:  graph.hpp
 *
 *    Description:  Class Graph header file
 *
 * 		  Version:  1.0
 *		  Created:  06/09/2016 16:22:24
 *		 Revision:  none
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


#include <unordered_set>
#include <map>

using namespace std;

typedef unordered_set<int> setOfEdges;

class Graph {

	private:
		int numOfVertex = 10;
		map<int, setOfEdges> G;

	public:
		Graph() {}
		Graph(int num): numOfVertex(num) {}
		
		map<int, setOfEdges>::iterator it;

		void addEdge(int src, int dest);
		void genRandomGraph(int num);
		void print();
};

#endif
