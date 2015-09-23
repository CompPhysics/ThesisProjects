#ifndef Variations_IS_INCLUDED
#define Variations_IS_INCLUDED

#include "../Domain/Domain.h"
#include "../Correlation/Correlation.h"
#include "../SlaterDet/SlaterDet.h"
#include "../Ref/Ref.h"
#include "../fFunction/fFunction.h"
#include "../LocalWaveFunction/LocalWaveFunction.h"
#include "../SpinFactors/SpinFactors.h"

// ****************************************************************
// *                          VARIATIONS                          *
// ****************************************************************
template <class SlaterDeterminant>
class Variations {

 protected:
  // !!! Temporary !!!
  ofstream  distanceFile;
  int       distanceFileIndex;


  ofstream* outputFile;
  fFunction*  f;

  int       numSlaterVariations;   // Number of variations in the
                                   // Slater-determinant parameter 
                                   // subspace. 
  int       numCorrelVariations;   // Number of variations in the
                                   // correlation parameter subspace.
  int       numVariations;         // Total number of variations.
  int       numParticles;          // Number of particles.

  int       centralSlater;         // Index specifying which of the 
                                   // Slater-variations we vary about.
  int       centralCorrel;         // Index specifying which of the 
                                   // correlation-variations we vary 
                                   // about.

  int       allowSpinFlip;         // Boolean indicating whether we 
                                   // allow(1) spin-flip or not(0).
  int       varianceOptimization;  // Boolean indicating whether we
                                   // optimize the paramters with 
                                   // respect to energy or variance
                                   // minimization.

  Ref<Domain>            domain;
  Ref<Correl>            correlation;
  Ref<SlaterDeterminant> slaterDeterminant;

  LocalWaveFunction*     localWaveFunctions;
  Correl*                correl;
  SlaterDeterminant*     slater;

  void                   createCorrel();
  void                   createSlaterDet();
  void                   createLocalSurface();

 public:
  Variations(Domain& _domain);
  void   setToNextParticle();
  void   setCurrentParticle(int _currentParticle);

  void   initializeThermalization(); // Called prior to starting the 
                                     // thermalization procedure.
  void   initializeVMC();            // Called in-between thermalization
                                     // and the actual VMC.
  void   initNewVmcRun();            // Called when subsequent VMC runs
                                     // are performed (after one VMC run
                                     // is finished and before next 
                                     // thermalization).

  void   suggestMove();
  void   acceptThermalization();
  void   rejectThermalization();
  void   acceptMove();
  void   rejectMove();

  void   sample();                   // One sample of the energy.
  void   summary();                  // Prints all the different local
                                     // variations to outputFile.
  void   summaryLowest();            // Prints all the lowest local
                                     // variation to outputFile (with
                                     // lowest we mean either lowest 
                                     // energy or variance).

  SlaterDeterminant& getSlaterDet()        {return slaterDeterminant();}
  Correl&            getCorrelation()      {return correlation();}

  int    findLowestWaveFunction();   // Finds the index corresponding 
                                     // to the lowest local variation.
  double getLowestEnergy();          // Find the energy corresponding 
                                     // to the lowest local variation.
  double getLowestVariance();        // Find the variance corresponding 
                                     // to the lowest local variation.
  double *getLowestAlphaParams();    // Find the Slater-paramters
                                     // corresponding to the lowest 
                                     // local variation.
  double *getLowestBetaParams();     // Find the correlation-paramters
                                     // corresponding to the lowest 
                                     // local variation.
  void setReferenceEnergy(double E); // Changes reference energy to E.

};

#include "Variations.cpp"


#endif
