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
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
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
    int cut_edge_num = 0;
    if (g.subgraphsNum() == 1) {
        return 0.0;
    }
    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        int vertex_colour = g.getColour(it->first);
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

void Analysis::cutEdgeVertexTable(const Graph& g, const vector<double>& ritz_values) {
    int subgraphs = g.subgraphsNum();

    std::vector<std::vector<int>> cut_edge_table(subgraphs, std::vector<int>(subgraphs, 0));
    std::vector<int> cut_vertex_table(subgraphs, 0);
    std::vector<int> isolated_vertex;

    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        int temp = 0;
        int vertex_subgraph = g.getColour(it->first);
        if (vertex_subgraph >= subgraphs) {
            vertex_subgraph = (vertex_subgraph/subgraphs)%subgraphs;
        }
        for (const int& neighbour:it->second) {
            int neighbour_subgraph = g.getColour(neighbour);
            if (neighbour_subgraph >= subgraphs) {
                neighbour_subgraph = (neighbour_subgraph/subgraphs)%subgraphs;
            }
            cut_edge_table[vertex_subgraph][neighbour_subgraph]++;
            if (vertex_subgraph == neighbour_subgraph) {
                temp++;
            }
        }
        if (temp == 0) {
            isolated_vertex.push_back(it->first);
        }
        cut_vertex_table[vertex_subgraph]++;
    }
    std::ostream_iterator<double> it_double(std::cout, "\t");
    std::ostream_iterator<int> it_int(std::cout, "\t");

    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << " * Basic info of the graph" << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << "Vertices:  " << g.size() << endl;
    cout << "Edges:     " << g.edgesNum() << endl;
    cout << "Colours:   " << subgraphs << endl;
    cout << "Used Ritz values: ";
    copy(ritz_values.cbegin(), ritz_values.cend(), it_double);
    cout << endl << "Cut Edge Percent: " << cutEdgePercent(g) * 100 << "%" << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << " * Number of vertices in each subgraph" << endl;
    cout << "/*-----------------------------------------------------------------------------" << endl;
    cout << "Colour:    " << "\t";
    for (int col = 0; col < subgraphs; col++) {
        cout << col << "\t";
    }
    cout << endl;
    cout << "Vertices:  " << "\t";
    copy(cut_vertex_table.cbegin(), cut_vertex_table.cend(), it_int);
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
    //if (isolated_vertex.size() != 0) {
    //    cout << "VertexIndex(Colour)" << "\t" << "NeighboursSize" << endl;
    //    for (const int& vertex:isolated_vertex) {
    //        cout << vertex << "(" << g.getColour(vertex) << ")" << "\t\t\t" << g.find(vertex)->second.size() << endl;
    //    }
    //}
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
    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        g.setColour(it->first, colour);
        num++;
        if (num == subgraph_size) {
            colour++;
            num = 0;
        }
    }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  outputTimes
 *  Description:  output the times into a file
 * =====================================================================================
 */

void Analysis::outputTimes(const int& procs, const int& size, const vector<double>& vec) {
    cout << "In P0, Lanczos takes " << vec[0] << "s" << endl;
    cout << "In P0, TQLI takes " << vec[1] << "s" << endl;
    cout << "In P0, Partition takes " << vec[2] << "s" << endl;

    string filename("./times/");
    filename += to_string(procs);
    filename += ".dat";
    ofstream Output(filename, ios::out | ios::binary | ios::app);
    //Output << "procs\tlanczos\ttqli\tpartition\t" << endl;
    Output << procs << "\t" << size << "\t";
    std::ostream_iterator<double> outIter(Output, "\t");
    std::copy(vec.cbegin(), vec.cend(), outIter);
    Output << endl;
    Output.close();
}
