/**
 * @file graph.h
 * @brief The interface of a graph
 * @author Ken Hu, xnchnhu@gmail.com
 */

/*
 * =====================================================================================
 *        Class:  Graph
 *  Description:  Class to create a new graph object
 * =====================================================================================
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

class Graph
{
private:
    typedef std::unordered_set<int> SetOfNeighbours;
    std::unordered_map<int, SetOfNeighbours> G;
    mutable std::unordered_map<int, int> Colour;

public:
    Graph() {}
    Graph(int n);  // Construct a random graph with n vertices

    void addEdge(int src, int dest);
    const int edgesNum() const;
    const int subgraphsNum() const;
    const int size() const;

    void outputDotFormat(const std::string& filename) const;
    void printLaplacianMat() const;
    void setColour(int vertex, int colour) const;
    const int getColour(int vertex) const;
    const int globalIndex(int& vertex) const;
    void readDotFormat(const std::string& filename);
    void readDotFormatWithColour(const std::string& filename);

    typedef std::unordered_map<int, std::unordered_set<int>>::const_iterator
        const_iterator;
    const const_iterator find(int vertex) const;
    const const_iterator cbegin() const;
    const const_iterator cend() const;
};

#endif
