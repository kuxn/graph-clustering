/**
 * @file test.h
 * @brief Header file for test.cc
 * @author Ken Hu, xnchnhu@gmail.com
 */

#ifndef TEST_H_
#define TEST_H_

class Tests {
    public:
        static bool testTqli();
        static bool testLanczos();
        static bool testPartition();
        static bool testReadGraph();
        static bool testReadGraphWithColour();
        static bool testRandomPartition();
        static bool testEvenPartition();
        static bool testCutEdgeVertexTable();
        static bool testReothogonalisation();
};

#endif

