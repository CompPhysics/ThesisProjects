#include "SpinFactors.h"

void SpinFactors::allocate(int _numParticles, int* _spinArray) {
  numParticles    = _numParticles;
  Nmatrix         = (numParticles*(numParticles-1))/2;
  spinArray       = _spinArray;
  newFactor       = new double[numParticles - 1];
  oldFactor       = new double[numParticles - 1];
  otherDifference = new double[numParticles - 1];
  matrix          = new double[Nmatrix];
  init();
}

void SpinFactors::init() {
  setCurrentParticle(0);
  setOtherParticle(0);
  setSpinFlip(0);
  setNewFactor(1);
}


// Call this routine if no spins are flipped 
// (OR the two spins are equal) 
void SpinFactors::calculateNoFlipFactors() {
  _spinArray  = spinArray;
  currentSpin = _spinArray + currentParticle;
  _newFactor  = newFactor  - 1;
  for (int i=0; i<currentParticle; i++)
    if ( (*_spinArray++)==(*currentSpin) ) (*++_newFactor)=0.5;
    else (*++_newFactor)=1;
  for (int i=currentParticle+1; i<numParticles; i++)
    if ( (*++_spinArray)==(*currentSpin) ) (*++_newFactor)=0.5;
    else (*++_newFactor)=1;
}

// Call this routine ONLY if the two spins differ
void SpinFactors::calculateFlipFactors() {
  _spinArray  = spinArray;
  currentSpin = spinArray + currentParticle;
  _newFactor  = newFactor - 1;
  _oldFactor  = oldFactor - 1;
  for (int i=0; i<currentParticle; i++)
    if ( (*_spinArray++)==(*currentSpin) ) {
      (*++_newFactor)=1;
      (*++_oldFactor)=0.5;
    }
    else {
      (*++_newFactor)=0.5;
      (*++_oldFactor)=1;
    }
  for (int i=currentParticle+1; i<numParticles; i++)
    if ( (*++_spinArray)==(*currentSpin) ) {
      (*++_newFactor)=1;
      (*++_oldFactor)=0.5;
    }
    else {
      (*++_newFactor)=0.5;
      (*++_oldFactor)=1;
    }
  *(newFactor+otherParticle-(otherParticle>currentParticle)) = 1;

  _spinArray       = spinArray;
  otherSpin        = spinArray + otherParticle;
  _otherDifference = otherDifference - 1;
  for (int i=0; i<otherParticle; i++)
    if ( (*_spinArray++)==(*otherSpin) ) (*++_otherDifference)=0.5;
    else (*++_otherDifference)=-0.5;
  for (int i=otherParticle+1; i<numParticles; i++)
    if ( (*_spinArray++)==(*otherSpin) ) (*++_otherDifference)=0.5;
    else (*++_otherDifference)=-0.5;
  *(otherDifference+currentParticle-(currentParticle>otherParticle)) = 0;
}


void SpinFactors::setToNextParticle() {
  currentParticle++;
  if (currentParticle == numParticles) setCurrentParticle(0);
}


void SpinFactors::calculateMatrix() {
  _matrix = matrix-1;
  for (int i=0; i<numParticles-1; i++)
    for (int j=i+1; j<numParticles; j++)
      if ( ( spinArray[i]==spinArray[j] ) & (quarterCusp) ) 
	(*++_matrix) = 0.5;
      else (*++_matrix) = 1;
}


void SpinFactors::setNewFactor(double value) {
  _newFactor = newFactor - 1;
  for (int i=0; i<numParticles-1; i++)
    (*++_newFactor) = value;
}


void SpinFactors::dumpNewFactor() {
  cerr << "Current Particle=" << currentParticle << " newFactor=[ ";
  _newFactor = newFactor - 1;
  for (int i=0; i<numParticles-1;i++)
    cerr << (*++_newFactor) << " ";
  cerr << " ]" << endl;
}


void SpinFactors::dumpOldFactor() {
  cerr << "Current Particle=" << currentParticle << " oldFactor=[ ";
  _oldFactor = oldFactor - 1;
  for (int i=0; i<numParticles-1;i++)
    cerr << (*++_oldFactor) << " ";
  cerr << " ]" << endl;
}


void SpinFactors::dumpOtherDifference() {
  cerr << "Current Particle=" << currentParticle << "Other Particle=" << otherParticle << " otherDifference=[ ";
  _otherDifference = otherDifference - 1;
  for (int i=0; i<numParticles-1;i++)
    cerr << (*++_otherDifference) << " ";
  cerr << " ]" << endl;
}


void SpinFactors::summaryNoFlip() {
  dumpNewFactor();
}


void SpinFactors::summaryFlip() {
  dumpNewFactor();
  dumpOldFactor();
  dumpOtherDifference();
}


void SpinFactors::summayMatrix() {
  cerr << "Matrix =" << endl;
  _matrix = matrix-1;
  for (int i=0; i<numParticles-1; i++) {
    for (int j=0; j<i+1; j++)
      cerr << 0 << " ";
    for (int j=i+1; j<numParticles; j++)
      cerr << (*++_matrix) << " ";
  cerr << endl;
  }
}

/*
void SpinFactors::() {
}
*/


