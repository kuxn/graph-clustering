// Test framework include
#include "wvtest.h"

// Include any "../*.h" header files which are required for your code
#include "../lanczos.h"
#include "../graph.h"
#include <iostream>
using namespace std;


WVTEST_MAIN("lanczos algorithm tests - dot product/norm")
{
	vector<double> v1(5, 1);
	vector<double> v2(5, 1);
	
	WVPASSEQ(dot(v1,v2), 5);

	vector<double> v3(9, 1);
	WVPASSEQ(norm(v3), 3);
	
}

WVTEST_MAIN("lanczos algorithm tests - Graph * Vector")
{
	Graph G;
	G.addEdge(0,1);
	G.addEdge(0,4);
	G.addEdge(1,2);
	G.addEdge(2,3);
	G.addEdge(2,4);
	G.addEdge(3,4);

	vector<double> v1(2, 1);
}


