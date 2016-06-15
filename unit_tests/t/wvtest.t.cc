// Test framework include
#include "wvtest.h"

// Include any "../*.h" header files which are required for your code
#include "../lanczos.h"
#include "../graph.h"
#include "../tqli.h"
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

int test_tqli() {
    int size = 4;
	map<pair<int,int>, double> eigenvec;
	for(int i = 0; i < size; i++) {
		eigenvec[make_pair(i,i)] = 1;
	}

	vector<double> diagonal, subdiagonal;
	diagonal = {6, 8, 11, 13};
	subdiagonal = {8, 2, 1, 0};

	tqli(diagonal, subdiagonal, size, eigenvec);
    vector<double> result = {-1.20826, 15.608, 13.2909, 10.3094};
    
    for (int i = 0; i < size; i++) {
        if (diagonal[i] - result[i] > 1e-5)
            return 0;
    }
    return 1;
}

WVTEST_MAIN("tqli algorithm tests - eigenvalues")
{
    WVPASS(test_tqli());
}




