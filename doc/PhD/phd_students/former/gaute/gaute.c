#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "complex.h"
#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-6
//#define vector
//#define free_vector
//typedef struct FCOMPLEX {double r,i;} fcomplex;
// definition of functions
fcomplex aa,bb,cc,z0,dz;
int kmax,kount;
float *xp,**yp,dxsav;
float **d,*x;

void bsstep(double  y[], double dydx[], int nv, double *xx, double htry,
	    double eps, double yscal[], double *hdid, double *hnext,
		void (*derivs)(double, double [], double []));
void hypdrv(double s, double yy[], double dyyds[]);
void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z,
	    fcomplex *series, fcomplex *deriv);
void odeint(double ystart[], int nvar, double x1, double x2,
	    double eps, double h1, double hmin, int *nok, int *nbad,
	    void (*derivs)(double, double [], double []),
	    void (*rkqs)(double [], double [], int, double *, double, double,
	    double [], double *, double *, void (*)(double, double [], double [])));
fcomplex hypgeo(fcomplex , fcomplex , fcomplex , fcomplex );



// begin main program

int main()
{
              fcomplex aa1,bb1,cc1,zz1,dzz1, number1;
	      //        double rr;
	      //  rr= 12.34345
	
// eksempel på kall
	    aa1=Complex(1.0,0.0);bb1=Complex(1.0,0.0);
	    cc1=Complex(1.0,0.0);dzz1=Complex(1.0,0.0);
	    zz1=Complex(1.0,0.0);
	    //  fprintf(stdout,"%13.3e \n",rr); 
	           number1 = hypgeo(aa1, bb1, cc1,zz1);

// deretter kan du skrive det heile ut vha f.eks printf e.l.

        return 0;
}  // end main program 


fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z)
{
	int nbad,nok;
        int kmax,kount;
        double *xp,**yp,dxsav;
	fcomplex ans,y[3];
        fcomplex aa,bb,cc,z0,dz;
	double *yy;

	kmax=0;
	if (z.r*z.r+z.i*z.i <= 0.25) {
		hypser(a,b,c,z,&ans,&y[2]);
		return ans;
	}
	else if (z.r < 0.0) z0=Complex(-0.5,0.0);
	else if (z.r <= 1.0) z0=Complex(0.5,0.0);
	else z0=Complex(0.0,z.i >= 0.0 ? 0.5 : -0.5);
	aa=a;
	bb=b;
	cc=c;
	dz=Csub(z,z0);
	hypser(aa,bb,cc,z0,&y[1],&y[2]);
	yy=vector(1,4);
	yy[1]=y[1].r;
	yy[2]=y[1].i;
	yy[3]=y[2].r;
	yy[4]=y[2].i;
	odeint(yy,4,0.0,1.0,EPS,0.1,0.0001,&nok,&nbad,hypdrv,bsstep);
	y[1]=Complex(yy[1],yy[2]);
	free_vector(yy,1,4);
	return y[1];
}











