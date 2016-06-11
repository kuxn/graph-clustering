/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Main function
 *
 * 		  Version:  1.0
 *		  Created:  06/09/2016 18:12:52
 *		 Revision:  none
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include "graph.hpp"


int main() {

	Graph G;

	G.addEdge(0,1);
	G.addEdge(0,4);
	G.addEdge(1,2);
	G.addEdge(2,3);
	G.addEdge(2,4);
	G.addEdge(3,4);
	
	//G.printDotFormat();

	G.printLaplacianMat();

	return 0;
}


