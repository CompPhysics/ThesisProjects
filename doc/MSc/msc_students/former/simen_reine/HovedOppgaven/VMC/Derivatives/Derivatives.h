#ifndef Derivatives_IS_INCLUDED
#define Derivatives_IS_INCLUDED

#include "../Jastrow/Jastrow.h"
#include "../Domain/Domain.h"


class Derivatives {

 protected:
  // Number of particles  
  int numParticles;
  int Nm1;
  int Nmatrix;
  // The particle currently (proposed) moved
  int currentParticle;
  Jastrow* jastrow;
  JastrowDiff* jastrowPlus;
  JastrowDiff* jastrowMinus;
  double twoh, hh;
  double* derivatives;
  double* derivativesNewColumn;
  double* secondDerivatives;
  double* secondDerivativesNewColumn;
  double* fPlusNewColumn;
  double* fMinusNewColumn;
  double* fNewColumn;
  double* _derivatives;
  double* _secondDerivatives;
  double* _derivativesNewColumn;
  double* _secondDerivativesNewColumn;
  void    nextParticle();
  void    setCurrentParticle(int _currentParticle);
  int*    spinArray;
  int*    _spinArrayI;
  int*    _spinArrayJ;
  double* cuspDerivatives;
  double* _cuspDerivatives;
  double* cuspSecondDerivatives;
  double* _cuspSecondDerivatives;
  

 public:
  Derivatives();
  void attach(Domain* domain, Jastrow* _jastrow, JastrowDiff* _jastrowPlus, 
	      JastrowDiff* _jastrowMinus);
  void initialize();
  void updateDerivatives();
  void acceptDerivatives();
  void rejectDerivatives();

  void getDColumn(int column, double* Column);
  void getD2Column(int column, double* Column);
  double* getDNewColumn()  {return derivativesNewColumn;}
  double* getD2NewColumn() {return secondDerivativesNewColumn;}
  double* getDMatrix()     {return derivatives;}
  double* getD2Matrix()    {return secondDerivatives;}
  double* getDCuspMatrix(); 
  double* getD2CuspMatrix();
  void summary();
};

#endif
