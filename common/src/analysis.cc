/**
 * @file analysis.cc
 * @brief Partitioning analysis
 * @author Ken Hu, xnchnhu@gmail.com
 */

#include "analysis.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include "partition.h"

using namespace std;

/**
 * @brief The percentage of edges have been cut by partitioning
 * @param g The graph to be analysed
 * @return The percentage of the cut edges
 */
double Analysis::cutEdgePercent(const Graph& g)
{
    int numOfCutEdges = 0;
    if (g.subgraphsNum() == 1) {
        return 0.0;
    }
    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        int vertexColour = g.getColour(it->first);
        for (const int& neighbour : it->second) {
            int neighbourColour = g.getColour(neighbour);
            if (neighbourColour != vertexColour) {
                numOfCutEdges++;
            }
        }
    }
    return (double)numOfCutEdges / (double)g.edgesNum() / 2.0;
}

/**
 * cutEdgeVertexTable
 * @param g graph to be analysed
 * @param ritzValues
 */
void Analysis::cutEdgeVertexTable(const Graph& g,
                                  const vector<double>& ritzValues)
{
    int subgraphs = g.subgraphsNum();

    std::vector<std::vector<int>> cut_edge_table(
        subgraphs, std::vector<int>(subgraphs, 0));
    std::vector<int> cut_vertex_table(subgraphs, 0);
    std::vector<int> isolated_vertex;

    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        int temp = 0;
        int vertex_subgraph = g.getColour(it->first);
        if (vertex_subgraph >= subgraphs) {
            vertex_subgraph = (vertex_subgraph / subgraphs) % subgraphs;
        }
        for (const int& neighbour : it->second) {
            int neighbour_subgraph = g.getColour(neighbour);
            if (neighbour_subgraph >= subgraphs) {
                neighbour_subgraph =
                    (neighbour_subgraph / subgraphs) % subgraphs;
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

    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    cout << " * Basic info of the graph" << endl;
    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    cout << "Vertices:  " << g.size() << endl;
    cout << "Edges:     " << g.edgesNum() << endl;
    cout << "Colours:   " << subgraphs << endl;
    cout << "Used Ritz values: ";
    copy(ritzValues.cbegin(), ritzValues.cend(), it_double);
    cout << endl
         << "Cut Edge Percent: " << cutEdgePercent(g) * 100 << "%" << endl;
    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    cout << " * Number of vertices in each subgraph" << endl;
    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    cout << "Colour:    "
         << "\t";
    for (int col = 0; col < subgraphs; col++) {
        cout << col << "\t";
    }
    cout << endl;
    cout << "Vertices:  "
         << "\t";
    copy(cut_vertex_table.cbegin(), cut_vertex_table.cend(), it_int);
    cout << endl;
    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    cout << " * Edges table after partitioning" << endl;
    cout << " * Each element represents number of edges (inside)/between "
            "subgraphs"
         << endl;
    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    for (int col = 0; col < subgraphs; col++) {
        cout << "\t" << col;
    }
    cout << endl;
    for (int row = 0; row < subgraphs; row++) {
        cout << row << "\t";
        for (int col = 0; col < subgraphs; col++) {
            if (row == col)
                cout << "(" << cut_edge_table[row][col] / 2 << ")"
                     << "\t";
            else
                cout << cut_edge_table[row][col] << "\t";
        }
        cout << endl;
    }
    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    cout << " * Vertices that have no neighbours in the same subgraph" << endl;
    cout << "/*----------------------------------------------------------------"
            "-------------"
         << endl;
    cout << "There are " << isolated_vertex.size() << " such vertices." << endl;
    // if (isolated_vertex.size() != 0) {
    //    cout << "VertexIndex(Colour)" << "\t" << "NeighboursSize" << endl;
    //    for (const int& vertex:isolated_vertex) {
    //        cout << vertex << "(" << g.getColour(vertex) << ")" << "\t\t\t" <<
    //        g.find(vertex)->second.size() << endl;
    //    }
    //}
}

/**
 * @brief Manually set equal number of vertices a colour and then compare the
 *        cutEdgePercent with Partitioning algorithm for parallel
 * @param g The graph to be analysed
 */
void Analysis::manuallyPartition(const Graph& g)
{
    int size = g.size();
    int subgraphSize = size / g.subgraphsNum();

    int num = 0, colour = 0;
    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        g.setColour(it->first, colour);
        num++;
        if (num == subgraphSize) {
            colour++;
            num = 0;
        }
    }
}

/**
 * @brief output the times into a file for parallel
 * @param procs Number of processes
 * @param numOfVertices
 * @param vec Vector contains the time elapsed in Lanczos, TQLI, Partition
 */

void Analysis::outputTimes(const int& procs, const int& numOfVertices,
                           const vector<double>& vec)
{
    cout << "In P0, Lanczos takes " << vec[0] << "s" << endl;
    cout << "In P0, TQLI takes " << vec[1] << "s" << endl;
    cout << "In P0, Partition takes " << vec[2] << "s" << endl;

    string filename("./times/");
    filename += to_string(procs);
    filename += ".dat";
    ofstream Output(filename, ios::out | ios::app);
    // Output << "procs\tlanczos\ttqli\tpartition\t" << endl;
    Output << procs << "\t" << numOfVertices << "\t";
    std::ostream_iterator<double> outIter(Output, "\t");
    std::copy(vec.cbegin(), vec.cend(), outIter);
    Output << endl;
    Output.close();
}

/**
 * @brief Analysis for serial
 */
void Analysis::randomPartition(const Graph& g, const int& colours)
{
    int colour;
    for (auto it = g.cbegin(); it != g.cend(); ++it) {
        colour = rand() % colours;
        g.setColour(it->first, colour);
    }
}

void Analysis::evenPartition(const Graph& g, const int& colours)
{
    int size = g.size();
    int subgraph_size = size / colours;

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

/**
 * @brief output the times into a file
 * @param numOfVertices
 * @param vec Vector contains the time elapsed in Lanczos, TQLI, Partition
 */
void Analysis::outputTimes(const int& numOfVertices, const vector<double>& vec)
{
    cout << "Lanczos takes " << vec[0] << "s" << endl;
    cout << "TQLI takes " << vec[1] << "s" << endl;
    cout << "Partition takes " << vec[2] << "s" << endl;
    string filename("./times/1.dat");
    ofstream Output(filename, ios::out | ios::app);
    // Output << "procs\tlanczos\ttqli\tpartition\t" << endl;
    Output << 1 << "\t" << numOfVertices << "\t";
    std::ostream_iterator<double> outIter(Output, "\t");
    std::copy(vec.cbegin(), vec.cend(), outIter);
    Output << endl;
    Output.close();
}

/**
 * @brief output time data for different number of vertices/colours
 * @param enableGramSchmidt
 */
/*
void Analysis::benchmarks(bool enableGramSchmidt)
{
    vector<std::string> filenames = {
        "./test/par_test_1024.dot",   "./test/par_test_2048.dot",
        "./test/par_test_5120.dot",   "./test/par_test_10240.dot",
        "./test/par_test_102400.dot",
    };
    // vector<int> subgraphs = {2, 4, 8, 16, 32, 64};
    vector<int> subgraphs = {8};
    vector<std::vector<std::vector<double>>> output_times;
    cout << "filenames.size() = " << filenames.size() << endl;

    // output_times[filename][subgraphs][times]
    for (unsigned int i = 0; i < filenames.size(); i++) {
        Graph g;
        g.readDotFormat(filenames[i]);
        vector<std::vector<double>> subgraph_times;
        for (unsigned int j = 0; j < subgraphs.size(); j++) {
            Partition partition(g, subgraphs[j], enableGramSchmidt);
            outputTimes(g.size(), partition.times);
            subgraph_times.push_back(partition.times);
        }
        output_times.push_back(subgraph_times);
    }

    // Output times: Lanczos - 0; TQLI - 1; Partition - 2;
    for (int i = 0; i < 3; i++) {
        string filename("./times/serial_");
        filename += to_string(i);
        filename += ".dat";
        ofstream Output(filename, ios::out | ios::binary | ios::trunc);
        Output << "subgraphs\t100\t200\t500\t1k\t5k\t10k" << endl;
        for (unsigned int row = 0; row < subgraphs.size(); row++) {
            Output << subgraphs[row] << "\t";
            for (unsigned int col = 0; col < filenames.size(); col++) {
                Output << output_times[col][row][i] << "\t";
            }
            Output << endl;
        }
        Output.close();
    }
}
*/
