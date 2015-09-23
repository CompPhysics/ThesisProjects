#ifndef Vmc_IS_INCLUDED
#define Vmc_IS_INCLUDED

#include "../Domain/Domain.h"
#include "../Correlation/Correlation.h"
#include "../SlaterDet/SlaterDet.h"
#include "../Ref/Ref.h"
#include "../Variations/Variations.h"
#include "../Walker/Walker.h"

// ****************************************************************
// *                             VMC                              *
// ****************************************************************
template <class SlaterDeterminant>
class Vmc {
 protected:

  Ref<Domain>                          domain;
  Ref<SlaterDeterminant>               slaterDeterminant;
  Ref<Correl>                          correlation;
  Ref<Walker<SlaterDeterminant> >      walker;
  Ref<Variations<SlaterDeterminant> >  wf;
  Ref<Random2>                         randomMetro;

  int      numThermalization;  // Number of thermalization steps
  int      numCycles;          // Number of Monte Carlo cycles

  ofstream *outputFile;        // Outputfile
  string   thermalizationType; // Determines way to perform 
                               // the thermalization(s).
  string   vmcType;            // Determines way to perform 
                               // the VMC run(s).
  int      uniDirectionalMovement; // Boolean determining whether
                               // we to search for a minima until
                               // movement in parameter space is 
                               // no longer uni-directional(1) or
                               // not(0).
  int      numberVmcRuns;      // Either the number of VMC runs 
                               // performed or the maximum number 
                               // of VMC runs we allow when 
                               // performing a (uni-directional) 
                               // search for minima.

  int      numberOfUniDirectionalMoves; // Number of uni-directional
                               // moves.
  int      rank;               // The rank or prosess number of the
                               // MPI run.
  double   centerRank;         // (Number of ranks + 1)/2.

 public:
  Vmc(Domain& _domain);
  void thermalization();
  void adaptiveStepThermalization(); 
  void vmcSome(); 
  void vmcOneParticleAtATime(); 
  void initNextRun();
  void run();
  
};

#include "Vmc.cpp"

#endif
