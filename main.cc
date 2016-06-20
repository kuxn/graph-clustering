/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  Main function
 *		  Created:  06/09/2016 18:12:52
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"

using namespace std;

int main() {

	Graph g;

	// Addedge function test
/*	
	g.addEdge(0,1);
	g.addEdge(0,4);
	g.addEdge(1,2);
	g.addEdge(2,3);
	g.addEdge(2,4);
	g.addEdge(3,4);
*/	
/*
	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(1,2);
	g.addEdge(2,3);
	g.addEdge(3,4);
	g.addEdge(3,5);
	g.addEdge(3,6);
	g.addEdge(3,7);
	g.addEdge(5,6);
	g.addEdge(4,6);
	g.addEdge(4,7);
*/
	g.genRandomGraph(6);
	//g.printDotFormat();
    g.printLaplacianMat();
  
    multiPartition(g);

	return 0;
}


