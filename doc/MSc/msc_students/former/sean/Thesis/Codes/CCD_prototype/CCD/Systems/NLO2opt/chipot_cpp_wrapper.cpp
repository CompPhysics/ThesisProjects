#include <iostream>
#include <stdio.h>
#include "chipot_cpp_wrapper.hpp"
#include <complex>
#include "SPBasis.hpp"

using namespace std;

extern "C" {
	// get real and imaginary components of matrix element.
	void chipot_f90_wrapper_(double *matel_real, double *matel_im,
							int *Nparticles, double *rho,
							int *ps, int *pt, int *ppx, int *ppy, int* ppz,
							int *qs, int *qt, int *qpx, int *qpy, int *qpz,
							int *rs, int *rt, int *rpx, int *rpy, int *rpz,
							int *ss, int *st, int *spx, int *spy, int *spz);
}

complex<double> chipot_cpp_wrapper(SPBasis *basis, double density, int p, int q, int r, int s){
// I think tz = -1 is neutrons for me, but +1 in the fortran code.
	printf("Checking states: p: %d, q: %d, r: %d, s: %d\n",p,q,r,s);
	
	int ps, pt, ppx, ppy, ppz;
	int qs, qt, qpx, qpy, qpz;
	int rs, rt, rpx, rpy, rpz;
	int ss, st, spx, spy, spz;

	ppx = basis->indexMap[p][0];
	ppy = basis->indexMap[p][1];
	ppz = basis->indexMap[p][2];
	ps = basis->indexMap[p][3];
	pt = 1*basis->indexMap[p][4]; // Switch sign for fortran code
	qpx = basis->indexMap[q][0];
	qpy = basis->indexMap[q][1];
	qpz = basis->indexMap[q][2];
	qs = basis->indexMap[q][3];
	qt = 1*basis->indexMap[q][4]; // Switch sign for fortran code
	rpx = basis->indexMap[r][0];
	rpy = basis->indexMap[r][1];
	rpz = basis->indexMap[r][2];
	rs = basis->indexMap[r][3];
	rt = 1*basis->indexMap[r][4]; // Switch sign for fortran code
	spx = basis->indexMap[s][0];
	spy = basis->indexMap[s][1];
	spz = basis->indexMap[s][2];
	ss = basis->indexMap[s][3];
	st = 1*basis->indexMap[s][4]; // Switch sign for fortran code

	int Nparticles = basis->nParticles;
	double rho = density;
	double matel_real, matel_im;

	printf("Printing quantum numbers for state: %d\n",p);
	printf("%d %d %d %d %d\n",ppx,ppy,ppz,ps,pt);
	printf("Printing quantum numbers for state: %d\n",q);
	printf("%d %d %d %d %d\n",qpx,qpy,qpz,qs,qt);
	printf("Printing quantum numbers for state: %d\n",r);
	printf("%d %d %d %d %d\n",rpx,rpy,rpz,rs,rt);
	printf("Printing quantum numbers for state: %d\n",s);
	printf("%d %d %d %d %d\n",spx,spy,spz,ss,st);
	printf("Nparticles: %d, rho: %.15f\n",Nparticles,rho);

	complex<double> result;
	chipot_f90_wrapper_(&matel_real, &matel_im,
						&Nparticles, &rho,
						&ps, &pt, &ppx, &ppy, &ppz,
						&qs, &qt, &qpx, &qpy, &qpz,
						&rs, &rt, &rpx, &rpy, &rpz,
						&ss, &st, &spx, &spy, &spz);
	result.real(matel_real);
	result.imag(matel_im);
	return result;
}

void chipot_regression_test(){
	int ps = 1;
	int pt = 1;
	int ppx = 0;
	int ppy = 0;
	int ppz = 0;
	int qs = -1;
	int qt = 1;
	int qpx = 1;
	int qpy = 0;
	int qpz = 0;
	int rs = 1;
	int rt = 1;
	int rpx = 1;
	int rpy = 0;
	int rpz = 0;
	int ss = -1;
	int st = 1;
	int spx = 0;
	int spy = 0;
	int spz = 0;
	int Nparticles = 4;
		//double rho = 0.08;
	// I'm going to set rho equal to the fortran value to keep things consistent for now.
	double rho = 7.9999998211860657E-002; 
	// Wow, floating point is crazy, when I set rho = 0.08 in the fortran code and print it, it prints as 7.9999998211860657E-002, but when I set it to 0.08 in this wrapper, and print it in the fortran code, it prints as 8.0000000000000002E-002 which yields a matel_real: -25.466340086077199, which differs at the 7th sigfig! This makes sense since the fortran 0.08 fails after 7 sig figs. This is kinda fucked, since in the fortran code, rho is a real(wp) with wp = selected_real_kind(15). I would assume this would guarantee 0.08 stays as 0.08 until at least 15 sig figs, but I guess not...
	double matel_real = 0.0;
	double matel_im = 0.0;

	chipot_f90_wrapper_(&matel_real, &matel_im,
						&Nparticles, &rho,
						&ps, &pt, &ppx, &ppy, &ppz,
						&qs, &qt, &qpx, &qpy, &qpz,
						&rs, &rt, &rpx, &rpy, &rpz,
						&ss, &st, &spx, &spy, &spz);
	printf("In original fortran code:\nmatel:  ( -25.466339388117557     ,  0.0000000000000000     )\n");
	printf("Now: matel: (%.15f,%.15f)\n",matel_real,matel_im);

}
