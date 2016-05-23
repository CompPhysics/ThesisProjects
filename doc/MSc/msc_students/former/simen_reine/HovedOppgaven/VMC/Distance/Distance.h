#ifndef Distance_IS_INCLUDED
#define Distance_IS_INCLUDED

#include "../Coor/Coor.h"
#include <blitz/array.h>
#include <cmath>
using namespace blitz;

#ifndef sqr_IS_INCLUDED
#define sqr_IS_INCLUDED
inline double sqr(double x) {return x*x;};
#endif

// ****************************************************************
// *                          DISTANCE                            *
// ****************************************************************
class Distance {

 protected:
  CoorSpinDiff* Coordinate;
  CoorSpinDiff* TrialCoordinate;
  int numParticles;     // Number of particles
  int Nm1;              // numParticles - 1
  int numDimensions;    // Number of dimension
  // Array of Upper matrix composed of the distance between 
  // particle i and j
  Array<double, 2> interElectronicDistances;
  // The (new) distance from the current particle to the other particles
  Array<double, 1> trialColumn;
  Array<double, 1> trialRow;
  int currentParticle;  // This is the particle that currently is 
                        // (proposed) moved

 public:
  Distance() {}
  void    attach(CoorSpinDiff* _Coordinate, CoorSpinDiff* _TrialCoordinate, 
		 int _numParticles);
  void    initialize();
  void    setToNextParticle();
  void    setCurrentParticle(int _currentParticle);
  void    suggestMove();
  void    acceptMove();
  void    rejectMove() {}
  double  getPotential();
  Array<double, 2> getInterElectronicDistances() 
    {return interElectronicDistances;}
  Array<double, 1> getTrialColumn() {return trialColumn;}
  Array<double, 1> getTrialRow()    {return trialRow;}
};

// ****************************************************************
// *                        DISTANCEDIFF                          *
// ****************************************************************
class DistanceDiff : public Distance {
 protected:
  double h;          // Differential parameter, 
                     // dr/dx ~= ( r(x+h) - r(x-h) )/2h
  double twoh;       // 2*h
  int differentiate; // Which dimension (x=0, y=1 or z=2 in three dim.) 
                     // to be differentiated


 public:
  DistanceDiff() {}
  void attach(CoorSpinDiff* _Coordinate, CoorSpinDiff* _TrialCoordinate, 
	      int _numParticles, double _h, int _differentiate);
  void initialize();
  void suggestMove();
  void acceptMove();
  void rejectMove() {}
};

#endif
