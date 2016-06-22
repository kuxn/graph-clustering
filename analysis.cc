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

#include "analysis.h"
#include "partition.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cutEdgePercent
 *  Description:  The percentage of edges have been cut by partitioning
 * =====================================================================================
 */

double cutEdgePercent(const Graph& g) {
	//partition(g);

	int size = g.size();
	int cut_edge_num = 0;
	int total_edge_num = 0;

	for (int vertex = 0; vertex < size; vertex++) {
		int vertex_colour = g.getColour(vertex);
		auto it = g.find(vertex);
		total_edge_num += it->second.size();

		for (const int& neighbour:it->second) {
			int neighbour_colour = g.getColour(neighbour);
			if (neighbour_colour != vertex_colour) {
				cut_edge_num++;
			}
		}	
	}
	return (double)cut_edge_num/(double)total_edge_num; 
}



