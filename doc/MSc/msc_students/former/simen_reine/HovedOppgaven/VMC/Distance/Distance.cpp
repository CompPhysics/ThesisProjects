#include "Distance.h"

// ****************************************************************
// *                          DISTANCE                            *
// ****************************************************************
//
//
// Class to monitor inter-electronic distances:
//
// Here the interElectronicDistances matrix is given by:
//                    
//                  |  0  r01 r02 ... r0(N-1)  |
//                  | r01  0  r12 ... r1(N-1)  |
//                  | r02 r12  0  ... r2(N-1)  |
//                  |  .         .      .      |
//                  |  .           .    .      |
//                  |  .             .  .      |
//                  |  .           r(N-2)(N-1) |
//                  | r0(N-1) .  .  .   0      |
//                    
void Distance::attach(CoorSpinDiff* _Coordinate, 
		      CoorSpinDiff* _TrialCoordinate, 
		      int _numParticles) {
  Coordinate      = _Coordinate; 
  TrialCoordinate = _TrialCoordinate;
  numParticles    = _numParticles; 
  Nm1             = numParticles - 1;
  numDimensions   = (*Coordinate).getLen();

  interElectronicDistances.resize(numParticles, numParticles);
  trialColumn.resize(numParticles);
  trialRow.resize(numParticles);
}

void Distance::initialize() {
  CoorSpinDiff* temp;
  temp = TrialCoordinate;
  setCurrentParticle(0);
  for (int i=0; i<numParticles; i++) {
    TrialCoordinate = &Coordinate[i];
    suggestMove();
    acceptMove();
    setToNextParticle();
  }
  TrialCoordinate = temp;
}


void Distance::setToNextParticle() {
 currentParticle++;
 if (currentParticle == numParticles) currentParticle=0;
}


void Distance::setCurrentParticle(int _currentParticle) {
  currentParticle = _currentParticle;
}


void Distance::suggestMove()
{
  for (int i=0; i<numParticles; i++) {
    double difference = sqr(Coordinate[i]()-(*TrialCoordinate)());
    for (int j=1; j<numDimensions; j++) {
      Coordinate[i]++; (*TrialCoordinate)++;
      difference += sqr(Coordinate[i]()-(*TrialCoordinate)());
    }
    Coordinate[i].resetPtr(); (*TrialCoordinate).resetPtr();
    trialColumn(i) = trialRow(i) = sqrt(difference);
  }
  trialColumn(currentParticle) = trialRow(currentParticle) = 0;
  
}

void Distance::acceptMove()
{
  Range N(0,Nm1);
  interElectronicDistances( currentParticle, N )  = trialColumn( N );
  interElectronicDistances( N , currentParticle ) = trialRow( N );
}

double Distance::getPotential()
{
  double potential = 0;
  for (int i=0; i<Nm1; i++)
    for (int j=i+1; j<numParticles; j++) 
      potential += 1/interElectronicDistances(i, j);
  return potential;
}



// ****************************************************************
// *                        DISTANCEDIFF                          *
// ****************************************************************
//
//
// Class to monitor the values used for numerical derivation of the 
// Jastrow-factor:
//                    dfij/dxi     ~= ( fij(xi+h) - fij(xi-h))/2h
// and:
//                    d^2fij/dxi^2 ~= ( fij(xi+h) + fij(xi-h) -2fij)/h^2
//
// ie. we need the values rij(xi+h) and rij(xi-h). Here xi is the either 
// x,y,z (for the three dimensional problem) of particle i.
// The relation:
//                    drij/dxi = - drij/dxj
// Implies:
//                    rij(xi+h) = rij(xj-h)
//
// Utilizing this relation, we need only to calculate one of the two.
// Here the rij matrix, given by:
//                    
//        | 0         r01(x1+h) r02(x2+h) ... r0(N-1)(x(N-1)+h)   |
//        | r01(x0+h)    0      r12(x2+h) ... r1(N-1)(x(N-1)+h)   |
//        | r02(x0+h) r12(x2+h)     0     ... r2(N-1)(x(N-1)+h)   |
//        | .            .                            .           |
//        | .            .                            .           |
//        | .            .                            .           |
//        | .            .                  r(N-2)(N-1)(x(N-1)+h) |
//        | r0(N-1)(x0+h)          r(N-2)(N-1)(x(N-2)+h)    0     |
void DistanceDiff::attach(CoorSpinDiff* _Coordinate, 
			  CoorSpinDiff* _TrialCoordinate, 
			  int _numParticles, double _h, 
			  int _differentiate) {
  Coordinate      = _Coordinate; 
  TrialCoordinate = _TrialCoordinate;
  numParticles    = _numParticles; 
  Nm1             = numParticles - 1;

  numDimensions   = (*Coordinate).getLen();
  h               = _h;
  twoh            = 2*h;
  differentiate   = _differentiate;   

  interElectronicDistances.resize(numParticles, numParticles);
  trialColumn.resize(numParticles);
  trialRow.resize(numParticles);
}

void DistanceDiff::initialize() {
  CoorSpinDiff* temp;
  temp = TrialCoordinate;
  setCurrentParticle(0);
  for (int i=0; i<numParticles; i++) {
    TrialCoordinate = &Coordinate[i];
    suggestMove();
    acceptMove();
    setToNextParticle();
  }
  TrialCoordinate = temp;
}

void DistanceDiff::suggestMove()
{
  (*TrialCoordinate)(differentiate)-=h;
  for (int i=0; i<currentParticle; i++) {
    double difference = sqr(Coordinate[i]()-(*TrialCoordinate)());
    for (int j=1; j<numDimensions; j++) {
      Coordinate[i]++; (*TrialCoordinate)++;
      difference += sqr(Coordinate[i]()-(*TrialCoordinate)());
    }
    Coordinate[i].resetPtr(); (*TrialCoordinate).resetPtr();
    trialRow(i) = sqrt(difference);
  }
  for (int i=currentParticle+1; i<numParticles; i++) {
    double difference = sqr(Coordinate[i]()-(*TrialCoordinate)());
    for (int j=1; j<numDimensions; j++) {
      Coordinate[i]++; (*TrialCoordinate)++;
      difference += sqr(Coordinate[i]()-(*TrialCoordinate)());
    }
    Coordinate[i].resetPtr(); (*TrialCoordinate).resetPtr();
    trialRow(i) = sqrt(difference);
  }

  (*TrialCoordinate)(differentiate)+=twoh;
  for (int i=0; i<currentParticle; i++) {
    double difference = sqr(Coordinate[i]()-(*TrialCoordinate)());
    for (int j=1; j<numDimensions; j++) {
      Coordinate[i]++; (*TrialCoordinate)++;
      difference += sqr(Coordinate[i]()-(*TrialCoordinate)());
    }
    Coordinate[i].resetPtr(); (*TrialCoordinate).resetPtr();
    trialColumn(i) = sqrt(difference);
  }
  for (int i=currentParticle+1; i<numParticles; i++) {
    double difference = sqr(Coordinate[i]()-(*TrialCoordinate)());
    for (int j=1; j<numDimensions; j++) {
      Coordinate[i]++; (*TrialCoordinate)++;
      difference += sqr(Coordinate[i]()-(*TrialCoordinate)());
    }
    Coordinate[i].resetPtr(); (*TrialCoordinate).resetPtr();
    trialColumn(i) = sqrt(difference);
  }
  (*TrialCoordinate)(differentiate)-=h;
  trialRow(currentParticle) = 0;
  trialColumn(currentParticle) = 0;
}

void DistanceDiff::acceptMove()
{
  Range N(0,Nm1);
  interElectronicDistances( currentParticle, N )  = trialColumn( N );
  interElectronicDistances( N, currentParticle ) = trialRow( N );
}
