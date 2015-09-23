#ifndef Correlation_IS_INCLUDED
#define Correlation_IS_INCLUDED

#include "../Domain/Domain.h"
#include "../Jastrow/Jastrow.h"
#include "../Distance/Distance.h"
#include "../fFunction/fFunction.h"
#include "../SpinFactors/SpinFactors.h"
#include <blitz/array.h>

// ****************************************************************
// *                           CORREL                             *
// ****************************************************************
class Correl {

 protected:

  Domain*       domain;
  int           numParticles;    // Number of particles.
  int           Nm1;             // numParticles - 1.
  int           numDimensions;   // Number of dimension.

  int           currentParticle; // This is the particle that 
                                 // currently is (proposed) moved.
  int           otherParticle;   // Particle to exchange spin with.

  Distance*     distance;
  DistanceDiff* distanceDiff;
  Jastrow*      jastrow;
  JastrowDiff*  jastrowDiff;

  int           expBool;         // 0 G=J, 1 G=exp(J)
  double*       beta;            // Variational parameters
  double        ratio;           // Gnew/Gold
  fFunction*    f;               // J=sum(i<j) f(r_ij, r_i, r_j) 
  double        h, hDoubleInverse, hSquaredInverse; // For numerical 
                                                    // derivation
  double*       gradJastrow;     // The gradient of G divided by G
  double        laplacian;       // The Laplacian of G divided by G
  double        correlation;     // The value of G

  SpinFactors*  spinFactors;
  double*       spinFactorMatrix;

  Array<double, 2> jastrowMatrix;
  Array<double, 2> *jastrowMatrixPlus, *jastrowMatrixMinus;


 public:
  Correl();
  void attach(Domain& _domain);
  void attachFFunction(fFunction* _f);
  void createJastrowAndJastrowDiff();

  void setToNextParticle();
  void setCurrentParticle(int _currentParticle); 
  void setOtherParticle(int _otherParticle) 
    {otherParticle=_otherParticle;}

  // Algorithms for proposing a Metropolis step
  void suggestMove();
  // Thermalization spesific algorithms; prior to samle
  void initializeThermalization();
  void acceptThermalizedMove();
  void rejectThermalizedMove();
  // VMC spesific algorithms
  void initializeVMC();
  void acceptMove();
  void rejectMove();
  // Only applicabel prior to accepting move; 
  // returns either Jnew/Jold or exp(Jnew)/exp(Jold)
  void calculateRatio();
  double* getRatioPtr()               {return &ratio;}

  // Only applicable after accepting move; 
  // returns either J or exp(J)
  void    calculateGradAndLaplacianRatios();
  double* getGradRatio()              {return gradJastrow;}
  double* getLaplaceRatio()           {return &laplacian;}
  void    calculateCorrelation();
  double* getCorrelationPtr()         {return &correlation;}
  double  operator()();
};


#endif
