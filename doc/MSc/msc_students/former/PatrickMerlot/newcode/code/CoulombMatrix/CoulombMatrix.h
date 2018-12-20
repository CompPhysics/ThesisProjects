#ifndef CoulombMatrix_h_IS_INCLUDED
#define CoulombMatrix_h_IS_INCLUDED

#include <iostream>
#include "singleParticleOrbitals.h"
#include <cmath>
#include <fstream>


// file        : CoulombMatrix.h
// description : definition of class CoulombMatrix

using namespace std;
  

class CoulombMatrix
{
private:
  const singleParticleOrbitals* Basis;// Pointer to a basis of single particle orbitals
  double ** Coulomb_Matrix;           // Coulomb Matrix of the system
  double lambda;                      // dimensionless constant characteristic of interaction strength

  void allocateMemory ();             // Allocate memory to the Coulomb matrix
  void fillingCoulombMatrix ();       // Computing the Coulomb Matrix elements
  
public:                // members visible also outside the class
  CoulombMatrix (const singleParticleOrbitals* S,const double lbd);// constructor
  ~CoulombMatrix ();                              // Destructor

 // given ch=2M+S, and the index of an element in this channel, read the index of the corresponding couple of states 
  int mappingCMelementToCoupleStates (const int channel, const int indexElement, const int CoupleStates);

 // given 4 states (i.e.their index), this function reads the corresponding Coulomb Matrix element
  double mappingSingleStatesToCMelement (const int indexState1, const int indexState2, const int indexState3, const int indexState4) const;
  
  // given the indices for 2 couple of states |i> and |j>, assuming i<j, read the corresponding Coulomb Matrix element
  double read_CoulombMatrixElement(const int channel, const int bra, const int ket) const;
  
  // Compute the coulomb matrix element <i|V|j>as
  double compute_CoulombMatrixElement (const int channel, const int i, const int j, FILE* f);
  
  // Compute the total number of element in the Coulomb Matrix
  int sizeCoulombMatrix ();

  // update the value of the interaction strength
  void setLambda(const double lbd);

  // Write to file the content of the sparse Coulomb matrix
  void save2file(ostream & mFile);
};

// compute a part of the Coulomb matrix element: <12|V|34>
double anisimovas (const int, const int, const int, const int, const int, const int, const int, const int);

//compute (-1)^k
int minusPower(const int);

// computes log(n!)
double LogFac (const int);

// Computes the first ratio in the Asinimovas expression
double LogRatio1(const int,const int,const int,const int);

// Computes the 2nd ratio in the Asinimovas expression
double LogRatio2(const int);

// computes first product of indices in the Anisimovas/Matulis expression
double Product1 (const int, const int, const int, const int, const int, const int, const int, const int);

// Computes the log of the 2nd product in the Asinimovas expression
double LogProduct2(const int,const int,const int,const int,const int,const int,const int,const int,const int,const int,const int,const int);

// Computes the log of the 3rd product in the Asinimovas expression
double LogProduct3(const int,const int,const int,const int,const int,const int,const int,const int);

// computes the factorial of a positive integer
int factorial(int);

// computes the factorial of a positive integer using log functions
double FastFactorial(int);

// The function gamma() computes the gamma function of a real x
double gamma(double);

// The function lgamma() computes the logarithm of the gamma function of real argument x
double lgamma(double x);


#endif

