// Test framework include
#include "wvtest.h"

// Include any "../*.h" header files which are required for your code
#include "../lanczos.h"
#include "../graph.h"
#include "../tqli.h"
#include "../test.h"
#include <iostream>

using namespace Tests;

using std::vector;
using std::unordered_map;
using std::map;

WVTEST_MAIN("lanczos algorithm tests - dot product/norm")
{
	vector<double> v1(5, 1);
	vector<double> v2(5, 1);
	
	WVPASSEQ(dot(v1,v2), 5);

	vector<double> v3(9, 1);
	WVPASSEQ(norm(v3), 3);
	
}

WVTEST_MAIN("tqli algorithm tests - eigenvalues")
{
    WVPASS(testTqli());
}


WVTEST_MAIN("lanczos algorithm tests - eigenvalues") 
{
	WVPASS(testLanczos());
}

