/*
 * =====================================================================================
 *
 *       Filename:  analysis.cc
 *
 *    Description:  Partitioning analysis
 *        Created:  06/22/2016 11:25:43
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include "analysis.h"
#include "partition.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cutEdgePercent
 *  Description:  The percentage of edges have been cut by partitioning
 * =====================================================================================
 */
using namespace std;

double cutEdgePercent(const Graph& g) {

	int subgraphs = 4;
	Partition partition(g, subgraphs);

	int size = g.size();
	int cut_edge_num = 0;

	for (int vertex = 0; vertex < size; vertex++) {
		int vertex_colour = g.getColour(vertex);
		auto it = g.find(vertex);

		for (const int& neighbour:it->second) {
			int neighbour_colour = g.getColour(neighbour);
			if (neighbour_colour != vertex_colour) {
				cut_edge_num++;
			}
		}	
	}
	return (double)cut_edge_num/(double)g.edgesNum()/2.0; 
}

/*-----------------------------------------------------------------------------
 *  Modified cutEdgePercent function calculating percentage of cut edges between different subgraphs
 *-----------------------------------------------------------------------------*/

void cutEdgeNodeTable(const Graph& g, const int& subgraphs) {
	int size = g.size();
	
	std::unordered_map<int, std::vector<int>> cut_edge_table;
	std::vector<int> cut_node_table(subgraphs, 0);
	std::vector<int> vinitial(subgraphs, 0);
	
	for (int i = 0; i < subgraphs; i++) {
		cut_edge_table[i] = vinitial;
	}
	
	for (int vertex = 0; vertex < size; vertex++) {
		int vertex_subgraph = g.getColour(vertex);
		auto it = g.find(vertex);
		for (const int& neighbour:it->second) {
			int neighbour_subgraph = g.getColour(neighbour);
			cut_edge_table[vertex_subgraph][neighbour_subgraph]++;
		}
		cut_node_table[vertex_subgraph]++;
	}	
	
	cout << "/*-----------------------------------------------------------------------------" << endl;
	cout << " * Number of nodes in each subgraph of a graph with " << g.size() << " nodes" << endl;
	cout << "/*-----------------------------------------------------------------------------" << endl;
	cout << "Colour:" << "\t";
	for (int col = 0; col < subgraphs; col++) {
		cout << col << "\t";
	}
	cout << endl;
	cout << "Nodes: " << "\t";
	for (const int& num:cut_node_table) {	
		cout << num << "\t";
	}
	cout << endl;
	cout << "/*-----------------------------------------------------------------------------" << endl;
	cout << " * Edges table after partitioning of a graph with " << g.edgesNum() << " edges" << endl;
	cout << " * Each element represents number of edges (inside)/between subgraphs" << endl;
	cout << "/*-----------------------------------------------------------------------------" << endl;
	for (int col = 0; col < subgraphs; col++) {
		cout << "\t" << col;
	}
	cout << endl;
	for (int row = 0; row < subgraphs; row++) {
		cout << row << "\t";
		for (int col = 0; col < subgraphs; col++) {
			if (row == col)
				cout << "(" << cut_edge_table[row][col]/2 << ")" << "\t";
			else
				cout << cut_edge_table[row][col] << "\t";
		}
		cout << endl;
    }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  manuallyPartition
 *  Description:  Manually set equal number of vertices a colour and then compare the cutEdgePercent with Partitioning algorithm
 * =====================================================================================
 */

void manuallyPartition(const Graph& g, const int& subgraphs) {
	int size = g.size();
	int subgraph_size = size/subgraphs;
	
	int num = 0, colour = 0;
	for (int vertex = 0; vertex < size; vertex++) {
		g.setColour(vertex, colour);
		num++;
		if (num == subgraph_size - 1) {
			colour++;
			num = 0;
		}
	}
}

