#ifndef fFunction_IS_INCLUDED
#define fFunction_IS_INCLUDED

#include <string>
#include <iostream>
#include <cstdio>
using namespace std;

// ****************************************************************
// *                          FFUNCTION                           *
// ****************************************************************
class fFunction {
protected:
  int     numberParams;  // Number of parameters
  double* params;        // Array of parameters
  int     expBool;       // If correlation G=J;           expBool=0
                         // Else if correlation G=exp(J); expBool=1

  virtual inline double   f(double rij, double ri, double rj)          
    {return 0;}
  virtual inline double   R(double r, int paramNumber) 
    {return 0;}
  double R_IJ, R_I, R_J;

public:
  fFunction() {}
  virtual void    attach(int _numberParams, double* _params, 
			 int _expBool) {
    numberParams = _numberParams;
    params       = _params;
    expBool      = _expBool;
  }
  virtual int     getExpBool()           {return expBool;}
  //  virtual double* getParams()            {return params;}
  //  virtual int     getNumberParams()      {return numberParams;}
  virtual double  operator()(double rij, double ri, double rj) 
    {return f(rij, ri, rj);}
};

//**************************** fNone ******************************
class fNone : public fFunction {
protected:
  virtual inline double f(double rij, double ri, double rj) {return 1.;}

public:
  fNone() {}
};

//**************************** fBeta ******************************
class fBeta: public fFunction {
protected:
  virtual inline double   f(double rij, double ri, double rj) 
    { return 0.5*rij/(1+params[0]*rij);}

public:
  fBeta() {}
};

//************************** fBetaMany ****************************
class fBetaMany: public fFunction {
protected:
  virtual inline double   f(double rij, double ri, double rj) 
    { return (0.5*rij + params[2]*rij*rij + params[4]*rij*rij*rij
	      + params[6]*rij*rij*rij*rij + params[8]*rij*rij*rij*rij*rij)
	/(1+params[0]*rij + params[1]*rij*rij + params[3]*rij*rij*rij 
	  + params[5]*rij*rij*rij*rij + params[7]*rij*rij*rij*rij*rij);}  

public:
  fBetaMany() {}
};

//**************************** fBeta2 *****************************
class fBeta2: public fFunction {
protected:
  virtual inline double   f(double rij, double ri, double rj) 
    {return params[1]*rij/(1+params[0]*rij);}
public:
  fBeta2() {}
};

//**************************** fBeta3 *****************************
class fBeta3: public fFunction {
protected:
  virtual inline double   f(double rij, double ri, double rj) 
    {return params[1]*rij/(1+params[0]*rij+params[2]*rij*rij);}
public:
  fBeta3() {}
};

//*************************** fBeta2r2 ****************************
class fBeta2r2: public fFunction {
protected:
  virtual inline double f(double rij, double ri, double rj) 
    {return rij/(1+params[0]*rij+params[1]*rij*rij);}

public:
  fBeta2r2() {}
};

//************************** fExtended ****************************
class fExtended: public fFunction {

protected:
  virtual inline double R(double r, int paramNumber) 
    {
      return r/(1+ params[paramNumber]*r);
    }

  virtual inline double f(double rij, double ri, double rj) 
    {
      R_IJ = R(rij, 0); R_I = R(ri, 1); R_J = R(rj, 1); 
      return 
	// e-e:                                 // m n o
	params[2]*R_IJ                          // 0 0 1
	//+ params[3]*R_IJ*R_IJ                   // 0 0 2

	// e-n:                                 // m n o
	+ params[4]*(R_I+R_J)                   // 1 0 0
	//+ params[5]*(R_I*R_I+R_J*R_J)           // 2 0 0

	// e-e-n:                               // m n o
	//+ params[6]*R_I*R_I*R_J*R_J             // 2 2 0
	//+ params[7]*R_IJ*R_IJ*(R_I*R_I+R_J*R_J) // 2 0 2
	;
    }

public:
  fExtended() {}
};

#endif
