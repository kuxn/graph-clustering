/*
 * =====================================================================================
 *
 *       Filename:  test.h
 *
 *    Description:  Header file for test.cc
 *        Created:   06/21/2016 16:38:27
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#ifndef TEST_H_
#define TEST_H_

class Tests {
    private:
        bool testReothogonalisation();
        bool testCutEdgeVertexTable();

    public:
        static bool testTqli();
        static bool testLanczos();
        static bool testPartition();
        static bool testReadGraph();
        static bool testReadGraphWithColour();
        static bool testManuallyPartition();
};

#endif

