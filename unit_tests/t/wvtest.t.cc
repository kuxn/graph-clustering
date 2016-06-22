// Test framework include
#include "wvtest.h"

// Include any "../*.h" header files which are required for your code
#include "../../include/lanczos.h"
#include "../../include/tqli.h"
#include "../../include/partition.h"
#include "../../include/test.h"

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
    WVPASS(Tests::testTqli());
}


WVTEST_MAIN("lanczos algorithm tests - eigenvalues") 
{
	WVPASS(Tests::testLanczos());
}

