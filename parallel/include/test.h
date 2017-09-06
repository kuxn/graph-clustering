/**
 * @file test.h
 * @brief Header file for test.cc
 * @author Ken Hu, xnchnhu@gmail.com
 */

#ifndef TEST_H_
#define TEST_H_

class Tests
{
public:
    static bool testPartition();
    static bool testReadGraph();
    static bool testReadByColour();
    static bool testPartitionWithClusters();
};

#endif
