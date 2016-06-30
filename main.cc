/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  Main function
 *        Created:  06/09/2016 18:12:52
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <iostream>
#include <sstream>
#include "graph.h"
#include "lanczos.h"
#include "tqli.h"
#include "partition.h"
#include "analysis.h"

using namespace std;

int main(int argc, char* argv[]) {

	/*-----------------------------------------------------------------------------
	 *  Test case for basic functions
	 *-----------------------------------------------------------------------------*/
	
	//g.addEdge(0,1);
	//g.addEdge(0,4);
	//g.addEdge(1,2);
	//g.addEdge(2,3);
	//g.addEdge(2,4);
	//g.addEdge(3,4);


	/*-----------------------------------------------------------------------------
	 *  Special test case for two partitioning
	 *-----------------------------------------------------------------------------*/
	
	// Have complex eigenvalue
	//g.addEdge(0,1);
	//g.addEdge(0,2);
	//g.addEdge(1,2);
	//g.addEdge(2,3);
	//g.addEdge(3,4);
	//g.addEdge(3,5);
	//g.addEdge(3,6);
	//g.addEdge(3,7);
	//g.addEdge(5,6);
	//g.addEdge(4,6);
	//g.addEdge(4,7);


	//g.addEdge(0,1);
	//g.addEdge(0,2);
	//g.addEdge(1,2);
	//g.addEdge(2,3);
	//g.addEdge(3,4);
	//g.addEdge(3,5);
	//g.addEdge(4,5);


	/*-----------------------------------------------------------------------------
	 *  Special case for multiple partition
	 *-----------------------------------------------------------------------------*/

	//g.addEdge(0,1);
	//g.addEdge(0,2);
	//g.addEdge(1,2);
	//g.addEdge(2,9);
	//g.addEdge(3,4);
	//g.addEdge(3,5);
	//g.addEdge(4,5);
	//g.addEdge(5,9);
	//g.addEdge(6,7);
	//g.addEdge(6,8);
	//g.addEdge(7,8);
	//g.addEdge(8,9);

	/*-----------------------------------------------------------------------------
	 *  Special case for multiple partition2
	 *-----------------------------------------------------------------------------*/

	//g.addEdge(0,1);
	//g.addEdge(0,2);
	//g.addEdge(1,2);
	//g.addEdge(2,12);
	//g.addEdge(3,4);
	//g.addEdge(3,5);
	//g.addEdge(4,5);
	//g.addEdge(5,12);
	//g.addEdge(6,7);
	//g.addEdge(6,8);
	//g.addEdge(7,8);
	//g.addEdge(8,12);
	//g.addEdge(9,10);
	//g.addEdge(9,11);
	//g.addEdge(10,11);
	//g.addEdge(10,12);
	//g.addEdge(13,12);
	//g.addEdge(13,14);
	//g.addEdge(13,15);
	//g.addEdge(14,15);

	istringstream ss(argv[1]);
	int num = 10;
	if (!(ss >> num))
		cerr << "Invalid number " << argv[1] << endl;

	cout << "num of vertices= " << num << endl;

	//Graph g(num);

	
	//g.printDotFormat();
	//g.outputDotFormat("test_8.dot");

	Graph g;
	ifstream In("par_test_200.dot");
	g.readDotFormatWithColour(In);
	//g.printLaplacianMat();
	Partition partition(g, 4, true);
	//partition.usingFullMat(g, 4, false);
	partition.printLapEigenvalues();
	//partition.printLapEigenMat();


	return 0;
}

