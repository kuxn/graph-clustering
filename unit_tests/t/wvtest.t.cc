// Test framework include
#include "wvtest.h"

// Include any "../*.h" header files which are required for your code
#include "../../include/lanczos.h"
#include "../../include/tqli.h"
#include "../../include/partition.h"
#include "../../include/test.h"
#include "../../include/analysis.h"

WVTEST_MAIN("tqli algorithm tests - eigenvalues")
{
    WVPASS(Tests::testTqli());
}

WVTEST_MAIN("lanczos algorithm tests - eigenvalues") 
{
	WVPASS(Tests::testLanczos());
}

WVTEST_MAIN("analysis - percentageo of cut edges") 
{
	WVPASS(Tests::testCutEdgePercent());
}


