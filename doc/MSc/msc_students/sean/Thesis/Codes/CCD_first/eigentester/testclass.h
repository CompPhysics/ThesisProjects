#ifndef TESTCLASS_H
#define TESTCLASS_H

#include <map>
#include <iostream>
#include <vector>

class testClass
{
public:
    testClass();

    void testFunc();
    int testy(int id);
    std::map<int, double> places;

    std::vector<int> D1;
    std::vector<int> D2;

};

#endif // TESTCLASS_H
