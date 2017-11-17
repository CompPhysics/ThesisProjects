#include <iostream>
#include <stdio.h>
#include "chipot_cpp_wrapper.hpp"
#include "SPBasis.hpp"
#include "InfMatterSPBasis.hpp"
#include <complex>

using namespace std;

SPBasis * basis;

//  Use makefile to compile
//  To run the regression test, use
//  ./main.exe 1 0.079999998211860657 2 2 1

int main(int argc, char * argv[]){

	int nParticles, nSpstates;
	double density;
	std::size_t basisIndicator, tzMax, nShells, nParticleShells;

	if( argc != 6){
		cout << "Please select : 0 xi g numPairStates nParticles: values at command line" << endl;
		cout << "or : 1 density tzMax shellMax nParticles: values at command line" << endl;
		cout << "or : 2 0 r_s shellMax nParticles: values at command line" << endl;
		cout << "or : 3 0 hbarOmega shellMax nParticles: values at command line" << endl;
		return 0;
	} else if ( atoi(argv[1]) == 1 ){
		basisIndicator = atoi(argv[1]);
		density = atof(argv[2]);
		tzMax = atof(argv[3]);
		nShells = atoi(argv[4]);
		nParticleShells = atoi(argv[5]);
		basis = new InfMatterSPBasis (basisIndicator,density, tzMax, nShells, nParticleShells);
	}

	nSpstates = basis->nSpstates;
	nParticles = basis->nParticles;
	//Check build
	//basis->printBasis();

	complex<double> answer;
	// For the regression test, check the matrix element
	// p = 3, q = 25, r = 27, s = 1
	answer = chipot_cpp_wrapper(basis, density, 3,25,27,1);

	printf("matel_real: %.15f\n",answer.real());
	printf("matel_im: %.15f\n",answer.imag());

	printf("Fortran Regression Test:\n");
	chipot_regression_test();
	return 0;
}
