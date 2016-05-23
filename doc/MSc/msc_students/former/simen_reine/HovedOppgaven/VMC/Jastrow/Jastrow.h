#ifndef Jastrow_IS_INCLUDED
#define Jastrow_IS_INCLUDED

#include "./../fFunction/fFunction.h"
#include "./../Distance/Distance.h"
#include "../SpinFactors/SpinFactors.h"
#include "../Domain/Domain.h"
#include "../Coor/Coor.h"
#include <blitz/array.h>

#include <iostream>
#include <cmath>
using namespace std;

class Jastrow {

 protected:
  Distance* distance;
  SpinFactors* spinFactors;

  fFunction* f;
  Domain* domain;

  CoorSpinDiff* coors;
  CoorSpinDiff* trialCoor;

  // Array of Upper matrix composed of the values of the function fij 
  Array<double, 2> jastrowMatrix;
  // The (new) values
  Array<double, 1> trialColumn, trialRow;
  Array<double, 1> trialDistanceColumn, trialDistanceRow;
  


  int numParticles;        // Number of particles.
  int Nm1;                 // numParticles - 1.
  int currentParticle;     // Particle proposed moved.
  double jastrowian;       // Value of J, G=J or G=exp(J).
  double difference;       // Value of deltaJ = Jnew - Jold.

  int* spinFlip;           // Boolean 0=noflip, 1=flip.
  int otherParticle;       // Particle to exchange spin with.
  double* newFactor;       // Values of either 1 or 1/2 to impose
                           // different cusp conditions.
  double* oldFactor;       // Same as above.
  double* otherDifference; // newFactor - oldFactor.

  void suggestMoveNoFlip();
  void calculateNoFlipDifference();
  void suggestMoveFlip();
  void calculateFlipDifference();

 public:
  Jastrow() {}
  void attach(Domain* _domain, fFunction* _f);
  void initialize();
  void setToNextParticle();
  void setCurrentParticle(int _currentParticle);
  void setOtherParticle(int _otherParticle) {otherParticle=_otherParticle;}
  void suggestMove();
  void acceptMove();
  void rejectMove() {}
  double& operator()();
  double getDifference();
  Array<double, 2> getJastrowMatrix() {return jastrowMatrix;}
  Array<double, 1> getTrialColumn()   {return trialColumn;}
  Array<double, 1> getTrialRow()      {return trialRow;}
};


class JastrowDiff{
  // The 2*d array, fij(xi_d+h) and fij(xi_d-h)
  DistanceDiff* distanceDiff;
  SpinFactors* spinFactors;
  fFunction* f;
  Domain* domain;
  CoorSpinDiff* coors;
  CoorSpinDiff* trialCoor;

  Array<double, 2> jastrowMatrixPlus, jastrowMatrixMinus;
  // The (new) values
  Array<double, 1> trialColumnPlus, trialRowPlus;
  Array<double, 1> trialColumnMinus, trialRowMinus; 
  Array<double, 1> trialDistanceColumn,  trialDistanceRow;

  int numParticles;    // Number of particles.
  int Nm1;             // numParticles - 1.
  int numDimensions;   // Number of dimensions.
  int currentParticle; // Particle proposed moved.

  int diffPlus;        // Index indicating with respect to which 
                       // Cartesian coordinate we wish to 
                       // perform positive differentiation
                       // (x(d=0), y(d=1), etc.) fij(xi_d+h).
                       // diffPlus = d.
  int diffMinus;       // Index indicating with respect to which 
                       // Cartesian coordinate we wish to
                       // perform nagative differentiation
                       // (x(d=0), y(d=1), etc.) fij(xi_d-h).
                       // diffPlus = d + numDimensions.

 public:
  JastrowDiff() {}
  void attach(Domain* _domain, fFunction* _f, int _differentiate);
  void initialize();
  void setToNextParticle();
  void setCurrentParticle(int _currentParticle);
  void suggestMove();
  void acceptMove();
  void rejectMove() {}
  Array<double, 2> getJastrowMatrixPlus()  {return jastrowMatrixPlus;}
  Array<double, 2> getJastrowMatrixMinus() {return jastrowMatrixMinus;}
};

#include "InlineJastrow.h"

#endif
