#include "TestCpp.h"
#include "MyArray.h"

// class definition
TestCpp::TestCpp(int N_, const MyArray<double>& A): N(N_){
	cout<<endl;
	cout << "Printing from the constructor of TestCpp" << endl;
	cout << "N: " << N << endl;

	someArray.redim(N, N);	// Now someArray = A and someArray.size = NxN

	fillMatrix(A);
}


void TestCpp::fillMatrix(const MyArray<double>& A){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			someArray(i,j) = A(i,j);
		}
	}
}


// print someArray
void TestCpp::print_(){
	
	cout<<endl;
	cout << "Printing from C++" << endl;

	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			cout << someArray(i,j) << endl;	
		}
	}
}
