#ifndef _TESTCPP_H_
#define _TESTCPP_H_

#include "MyArray.h"
#include <iostream>


using namespace std;

// class declaration
class TestCpp{

	private:
		int N;
		MyArray<double> someArray;	// initialize the default constructor
															  // someArray = NULL and someArray.size() = 0
		void fillMatrix(const MyArray<double>& A);

	public:
		TestCpp(int N_, const MyArray<double>& A);
		
		void print_();
};

#endif
