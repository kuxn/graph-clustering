// Test framework include
#include "wvtest.h"

// Include any "../*.h" header files which are required for your code
#include "../../serial/include/lanczos.h"
#include "../../serial/include/tqli.h"
#include "../../serial/include/partition.h"
#include "../../serial/include/test.h"
#include "../../serial/include/analysis.h"

WVTEST_MAIN("Graph properties - read dot format") {
    WVPASS(Tests::testReadGraph());
}

WVTEST_MAIN("Graph properties - read dot format with colour") {
    WVPASS(Tests::testReadGraphWithColour());
}

WVTEST_MAIN("Partition - TQLI") {
    WVPASS(Tests::testTqli());
}

WVTEST_MAIN("Partition - Lanczos & TQLI") {
    WVPASS(Tests::testLanczos());
}

WVTEST_MAIN("Partition - partition & cut edge percentage") {
    WVPASS(Tests::testPartition());
}

WVTEST_MAIN("Analysis - manually partition & cut edge percentage") {
    WVPASS(Tests::testRandomPartition());
}
