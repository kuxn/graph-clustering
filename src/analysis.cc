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

using namespace std;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cutEdgePercent
 *  Description:  The percentage of edges have been cut by partitioning
 * =====================================================================================
 */

double Analysis::cutEdgePercent(const Graph& g) {

    int size = g.size();
    int cut_edge_num = 0;

    if (g.subgraphsNum() == 1) {
        return 0.0;
    }

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

void Analysis::cutEdgeVertexTable(const Graph& g) {
    int size = g.size();
    int subgraphs = g.subgraphsNum();

    //cout << "subgraphs = " << subgraphs << endl;

    std::unordered_map<int, std::vector<int>> cut_edge_table;
    std::vector<int> cut_vertex_table(subgraphs, 0);
    std::vector<int> isolated_vertex;

    std::vector<int> vinitial(subgraphs, 0);
    for (int i = 0; i < subgraphs; i++) {
        cut_edge_table[i] = vinitial;
    }

    for (int vertex = 0; vertex < size; vertex++) {
        int temp = 0;
        int vertex_subgraph = g.getColour(vertex);
        if (vertex_subgraph >= subgraphs)
            vertex_subgraph = (vertex_subgraph/subgraphs)%subgraphs;
        auto it = g.find(vertex);
        for (const int& neighbour:it->second) {
            int neighbour_subgraph = g.getColour(neighbour);
            if (neighbour_subgraph >= subgraphs)
                neighbour_subgraph = (neighbour_subgraph/subgraphs)%subgraphs;
            cut_edge_table[vertex_subgraph][neighbour_subgraph]++;
            if (vertex_subgraph == neighbour_subgraph) {
                temp++;
            }
        }
        if (temp == 0) {
            isolated_vertex.push_back(vertex);
        }
        cut_vertex_table[vertex_subgraph]++;
    }	

    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << " * Basic info of the graph" << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << "Vertices:  " << g.size() << endl;
    cout << "Edges:     " << g.edgesNum() << endl;
    cout << "Subgraphs: " << subgraphs << endl;
    cout << "Cut Edge Percent: " << cutEdgePercent(g) << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << " * Number of vertices in each subgraph" << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << "Colour:    " << "\t";
    for (int col = 0; col < subgraphs; col++) {
        cout << col << "\t";
    }
    cout << endl;
    cout << "vertices:  " << "\t";
    for (const int& num:cut_vertex_table) {	
        cout << num << "\t";
    }
    cout << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << " * Edges table after partitioning" << endl;
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
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << " * Vertices that have no neighbours in the same subgraph" << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << "There are " << isolated_vertex.size() << " such vertices." << endl;
    if (isolated_vertex.size() != 0) {
        cout << "VertexIndex(Colour)" << "\t" << "NeighboursSize" << endl;
        for (const int& vertex:isolated_vertex) {
            cout << vertex << "(" << g.getColour(vertex) << ")" << "\t\t\t" << g.find(vertex)->second.size() << endl;
        }
    }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  manuallyPartition
 *  Description:  Manually set equal number of vertices a colour and then compare the cutEdgePercent with Partitioning algorithm
 * =====================================================================================
 */

void Analysis::manuallyPartition(const Graph& g) {
    int size = g.size();
    int subgraph_size = size/g.subgraphsNum();

    int num = 0, colour = 0;
    for (int vertex = 0; vertex < size; vertex++) {
        g.setColour(vertex, colour);
        num++;
        if (num == subgraph_size) {
            colour++;
            num = 0;
        }
    }
}