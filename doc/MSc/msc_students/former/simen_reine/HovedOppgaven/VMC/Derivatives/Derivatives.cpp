#include "Derivatives.h"

Derivatives::Derivatives() {
}

void Derivatives::attach(Domain* domain, Jastrow* _jastrow, 
			 JastrowDiff* _jastrowPlus, 
			 JastrowDiff* _jastrowMinus) {

  jastrow                    = _jastrow; 
  jastrowPlus                = _jastrowPlus;
  jastrowMinus               = _jastrowMinus;

  numParticles               = domain->getNumParticles();
  Nmatrix                    = (numParticles*(numParticles-1))/2;
  Nm1                        = numParticles - 1;

  derivatives                = new double[Nmatrix];
  secondDerivatives          = new double[Nmatrix];
  derivativesNewColumn       = new double[Nm1];
  secondDerivativesNewColumn = new double[Nm1];
 
 double h                   = domain->getH();
  twoh                       = 2*h;
  hh                         = h*h;

  spinArray                  = domain->getSpinArray();
  cuspDerivatives            = new double[Nmatrix];
  cuspSecondDerivatives      = new double[Nmatrix];
}

void Derivatives::initialize() {
  double* jastrowMatrix = (*jastrow).getFMatrix() - 1;
  double* jastrowPlusMatrix = (*jastrowPlus).getFMatrix() - 1;
  double* jastrowMinusMatrix = (*jastrowMinus).getFMatrix() - 1;

  _derivatives       = derivatives - 1;
  _secondDerivatives = secondDerivatives - 1;

  for (int i=0; i<Nmatrix; i++) {
    (*++_derivatives) = 
      ( (*++jastrowPlusMatrix) - (*++jastrowMinusMatrix) ) / twoh;
    (*++_secondDerivatives) = 
      ( (*jastrowPlusMatrix) + (*jastrowMinusMatrix) 
	- 2*(*++jastrowMatrix) ) / hh;
  }
  setCurrentParticle(0);
}


void Derivatives::updateDerivatives() {

  fNewColumn                  = (*jastrow).getFNewColumn() - 1;
  fPlusNewColumn              = (*jastrowPlus).getFNewColumn() - 1;
  fMinusNewColumn             = (*jastrowMinus).getFNewColumn() - 1;

  _derivativesNewColumn       = derivativesNewColumn - 1;
  _secondDerivativesNewColumn = secondDerivativesNewColumn - 1;
  for (int i=0; i<Nm1; i++) {
    (*++_derivativesNewColumn) = 
      ( (*++fPlusNewColumn) - (*++fMinusNewColumn) ) / twoh;
    (*++_secondDerivativesNewColumn) = 
      ( (*fPlusNewColumn) + (*fMinusNewColumn) - 2*(*++fNewColumn) ) / hh;
  }
}

void Derivatives::acceptDerivatives() {
  _derivativesNewColumn       = derivativesNewColumn - 1;
  _secondDerivativesNewColumn = secondDerivativesNewColumn - 1;
  _derivatives                = derivatives + currentParticle - 1;
  _secondDerivatives          = secondDerivatives + currentParticle - 1;
  int k=Nm1;
  for (int i=0; i<currentParticle; i++) {
    (*_derivatives) = (*++_derivativesNewColumn);
    (*_secondDerivatives) = (*++_secondDerivativesNewColumn);
    k--;
    _derivatives+=k;
    _secondDerivatives+=k;
  }
  for (int i=currentParticle+1; i<numParticles; i++) {
    (*++_derivatives) = (*++_derivativesNewColumn);
    (*++_secondDerivatives) = (*++_secondDerivativesNewColumn);    
  }
  nextParticle();
}

void Derivatives::rejectDerivatives() {
  nextParticle();
}

void Derivatives::nextParticle() {
  currentParticle++;
  if (currentParticle==numParticles) currentParticle=0;
}

void Derivatives::setCurrentParticle(int _currentParticle) {
  currentParticle = _currentParticle;
}

void Derivatives::getDColumn(int column, double* Column) {
  _derivatives = derivatives + column -1;
  double* _Column=Column-1;
  int k=Nm1;
  for (int i=0; i<column; i++) {
    (*++_Column) = (*_derivatives);
    k--;
    _derivatives += k;
  }
  for (int i=column+1; i<numParticles; i++)
    (*++_Column) = (*++_derivatives);
}

void Derivatives::getD2Column(int column, double* Column) {
  _secondDerivatives = secondDerivatives + column -1;
  double* _Column=Column-1;
  int k=Nm1;
  for (int i=0; i<column; i++) {
    (*++_Column) = (*_secondDerivatives);
    k--;
    _secondDerivatives += k;
  }
  for (int i=column+1; i<numParticles; i++)
    (*++_Column) = (*++_secondDerivatives);
}


double* Derivatives::getDCuspMatrix() {
  _cuspDerivatives = cuspDerivatives - 1;
  _derivatives     = derivatives - 1;
  for (int i=0; i<Nm1; i++) {
    _spinArrayJ = _spinArrayI = spinArray + i;
    for (int j=i+1; j<numParticles; j++)
      if ( (*++_spinArrayJ)==(*_spinArrayI) )
	(*++_cuspDerivatives) = 0.5*(*++_derivatives);
      else
	(*++_cuspDerivatives) = (*++_derivatives);
  }
  return cuspDerivatives;
}


double* Derivatives::getD2CuspMatrix() {

  return cuspSecondDerivatives;
}


void Derivatives::summary() {

  cerr << "\n\n--------------------------------------------------";
  cerr << "\n SUMMARY OF Derivatives INSTANCE\n\n";
  cerr << " Derivatives: " << endl;
  _derivatives = derivatives - 1;
  for (int i=0; i<(numParticles*(numParticles-1))/2; i++)
    cerr << (*++_derivatives) << " ";
  cerr << endl;
  _secondDerivatives = secondDerivatives - 1;
  for (int i=0; i<(numParticles*(numParticles-1))/2; i++)
    cerr << (*++_secondDerivatives) << " ";
  cerr << endl;


}




